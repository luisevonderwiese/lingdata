import lingdata.glottolog as glottolog
import lingdata.crawler as crawler
from lingdata.categorical import CategoricalData
import lingdata.pathbuilder as pb
import lingdata.params as params
import lingdata.native_data as native_data
import lingdata.state_encoding as state_encoding
import os
import pandas as pd
import traceback
from ete3 import Tree
import json
from ast import literal_eval

columns = [
            "ds_id",
            "source",
            "ling_type",
            "family",
            "sub_families",
            "num_taxa",
            "num_chars",
            "num_sites_bin"
            "sites_per_char",
            "multistate_ratio",
            "value_number_counts",
            "max_values",
            "MULTIx_MK",
            "MULTIx_GTR",
            "MULTIx_MK+M",
            "MULTIx_GTR+M",
            "categorical_path",
            "glottolog_tree_path",
            "msa_paths",
            "primary_source",
            "modification_date",
            "glottolog_version"
            ]

converters={"value_number_counts": lambda x: [int(el) for el in x.strip("[]").split(", ")],
                "sub_families": lambda x: [] if x == "set()" else [el.strip("'")  for el in x.strip("{}").split(", ")],
                "msa_paths" : literal_eval}

def read_config(config_path):
    with open(config_path, 'r') as openfile:
        json_object = json.load(openfile)
    params.size_limit = json_object["size_limit"]
    params.family_split_threshold = json_object["family_split_threshold"]

    params.data_dir = json_object["data_dir"]
    params.native_dir = json_object["native_dir"]

    params.sources = json_object["sources"]
    params.ling_types = json_object["ling_types"]

    params.msa_types = json_object["msa_types"]

    params.glottolog_tree_required = json_object["glottolog_tree_required"]

    params.flat_paths = json_object["flat_paths"]

def data():
    metadata_path = pb.metadata_path()
    if not os.path.isfile(metadata_path):
        return None
    else:
        return  pd.read_csv(metadata_path, sep = ";", converters = converters)



def build_database(data_units):
    db_df = pd.DataFrame(columns = columns)
    for i, data_unit in enumerate(data_units):
        #general
        db_df.at[i, "ds_id"] = data_unit.ds_id
        db_df.at[i, "source"] = data_unit.source
        db_df.at[i, "ling_type"] = data_unit.ling_type
        db_df.at[i, "family"] = data_unit.family
        #families
        db_df.at[i, "sub_families"] = data_unit.sub_families
        #categorical data properties
        data = data_unit.data
        db_df.at[i, "num_taxa"] = data.num_taxa()
        db_df.at[i, "num_chars"] = data.num_chars()
        db_df.at[i, "num_sites_bin"] = data.num_sites_bin()
        db_df.at[i, "sites_per_char"] = data.sites_per_char()
        db_df.at[i, "multistate_ratio"] = data.get_multistate_ratio()
        db_df.at[i, "value_number_counts"] = data.get_value_number_counts()
        x = data.max_values()
        db_df.at[i, "max_values"] = x
        #multi models
        db_df.at[i, "MULTIx_MK"] = "MULTI" + str(x) + "_MK"
        db_df.at[i, "MULTIx_GTR"] = "MULTI" + str(x) + "_GTR"
        if "ambig" in params.msa_types:
            charmap_path = pb.charmap_path(x)
            db_df.at[i, "MULTIx_MK+M"] = "MULTI" + str(x) + "_MK+M{" + charmap_path + "}"
            db_df.at[i, "MULTIx_GTR+M"] = "MULTI" + str(x) + "_GTR+M{" + charmap_path + "}"
        else:
            db_df.at[i, "MULTIx_MK+M"] = ""
            db_df.at[i, "MULTIx_GTR+M"] = ""
        #paths
        db_df.at[i, "categorical_path"] = data_unit.categorical_path()
        tree_path = data_unit.glottolog_tree_path()
        if os.path.isfile(tree_path):
            db_df.at[i, "glottolog_tree_path"] = tree_path
        else:
            db_df.at[i, "glottolog_tree_path"] = ""
        msa_paths = {}
        for msa_type in params.msa_types:
            msa_path = data_unit.msa_path(msa_type)
            if not os.path.isfile(msa_path):
                msa_path  = ""
            msa_paths[msa_type] = msa_path
        db_df.at[i, "msa_paths"] = msa_paths

        #source information
        db_df.at[i, "primary_source"] = data_unit.primary_source
        db_df.at[i, "modification_date"] = data_unit.date
        db_df.at[i, "glottolog_version"] = params.glottolog_version

    db_df.to_csv(pb.metadata_path(), sep = ";")





def generate_data():
    pb.clear_data()
    if not (os.path.isdir(pb.domain_path("glottolog")) and os.path.isdir(pb.domain_path("native"))):
        print("Native data does not exist, run update_native() before")
        return

    if "ambig" in params.msa_types:
        state_encoding.write_charmaps()
    data_units = []
    for ds_id in os.listdir(pb.domain_path("native")):
        print("Checking data for ", ds_id)
        for source in params.sources:
            for ling_type in params.ling_types:
                data_units += generate_data_units(ds_id, source, ling_type)
    build_database(data_units)


def generate_data_units(ds_id, source, ling_type):
    handler = native_data.get_handler(ds_id, source)
    data_list = handler.get_data(ling_type)
    data_units = []
    if data_list == []:
        return data_units

    #general
    date = handler.get_date()
    primary_source = handler.get_source()

    for data, family in data_list:
        if data is None:
            continue
        data_unit = DataUnit(data)
        data_unit.ds_id = ds_id
        data_unit.source = source
        data_unit.ling_type = ling_type
        data_unit.family = family
        data_unit.date = date
        data_unit.primary_source = primary_source


        #tree
        (glottocodes, complete) = handler.get_glottocodes(data.taxon_ids)
        tree_path = data_unit.glottolog_tree_path()
        if complete:
            glottolog_tree = glottolog.get_tree(glottocodes)
            if glottolog_tree is not None:
                pb.mk_file_dir(tree_path)
                glottolog_tree.write(outfile = tree_path, format=9)
            elif params.glottolog_tree_required:
                continue

        #families
        if family == "full":
            data_unit.sub_families = glottolog.get_families(glottocodes)
        else:
            data_unit.sub_families = [family]

        #categorical
        categorical_path = data_unit.categorical_path()
        pb.mk_file_dir(categorical_path)
        data.write(categorical_path)

        #msas
        pb.mk_this_dir(pb.family_path("msa", ds_id, source, ling_type, family))
        for msa_type in params.msa_types:
            data.write_msa(data_unit.msa_path(msa_type), msa_type)

        data_units.append(data_unit)
        print("Data created for ", ds_id, source, ling_type, family)

    return data_units

def update_native():
    glottolog.crawl()
    crawler.crawl()



class DataUnit:

    data = None

    ds_id = ""
    source = ""
    ling_type = ""
    family = ""

    primary_source = ""
    date = None
    sub_families = []


    def __init__(self, data):
        self.data = data

    def categorical_path(self):
        return pb.categorical_path(self.ds_id, self.source, self.ling_type, self.family)

    def msa_path(self, msa_type):
        return pb.msa_path(self.ds_id, self.source, self.ling_type, self.family, msa_type)

    def glottolog_tree_path(self):
        return pb.glottolog_tree_path(self.ds_id, self.source, self.family)
