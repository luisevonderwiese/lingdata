import os
import pandas as pd
import json
from ast import literal_eval
from getpass4 import getpass
from termcolor import colored
from datetime import datetime

import lingdata.glottolog as glottolog
import lingdata.crawler as crawler
import lingdata.pathbuilder as pb
import lingdata.params as params
import lingdata.native_data as native_data
import lingdata.state_encoding as state_encoding
import lingdata.membership as membership


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
            "site_group_sizes",
            "multistate_ratio",
            "value_number_counts",
            "value_number_matrix",
            "max_values",
            "max_values_prototype",
            "num_discarded_prototype",
            "MULTIx_MK",
            "MULTIx_GTR",
            "MULTIx_MK+M",
            "MULTIx_GTR+M",
            "MULTIx_MK_prototype",
            "MULTIx_GTR_prototype",
            "COGx",
            "categorical_path",
            "glottolog_tree_path",
            "msa_paths",
            "primary_source",
            "native_data_source_sha",
            "glottolog_version",
            "download_cutoff"
            ]

converters={"value_number_counts": lambda x: [int(el) for el in x.strip("[]").split(", ")],
            "site_group_sizes": lambda x: [int(el) for el in x.strip("[]").split(", ")],
           # "value_number_matrix": lambda x: [] if x == "[]" else [[] if el == "[]" else [int(inner_el) for inner_el in el.strip("[]").split(", ")] for el in x.strip("[]").split(", ")],
                "sub_families": lambda x: [] if x == "set()" else [el.strip("'")  for el in x.strip("{}").split(", ")],
                "sampled_msa_paths": lambda x: [] if x == "[]" else [el.strip("'") for el in x.strip("[]").split(", ")],
                "msa_paths" : literal_eval,
                "partition_paths" : literal_eval}

all_sources = ["lexibank", "SequenceComparison", "correspondence-pattern-data"]
all_ling_types = ["cognate", "structural", "correspondence"]
all_msa_types = ["bin", "multi", "catg_bin", "catg_multi", "ambig", "membership", "prototype"]

def read_config(config_path):
    with open(config_path, 'r') as openfile:
        json_object = json.load(openfile)
    params.max_num_taxa = json_object["max_num_taxa"]
    if not isinstance(params.max_num_taxa, int):
        raise Exception("Malformed config: max_num_taxa must be an integer")
    params.max_num_chars = json_object["max_num_chars"]
    if not isinstance(params.max_num_chars, int):
        raise Exception("Malformed config: max_num_chars must be an integer")
    params.family_split_threshold = json_object["family_split_threshold"]
    if not isinstance(params.family_split_threshold, int):
        raise Exception("Malformed config: family_split_threshold must be an integer")
    params.num_samples = json_object["num_samples"]
    if not isinstance(params.num_samples, int):
        raise Exception("Malformed config: num_samples must be an integer")

    config_dir = os.path.dirname(config_path)
    cur_cwd = os.getcwd()
    if config_dir != "":
        os.chdir(config_dir)
    params.data_dir = os.path.abspath(json_object["data_dir"])
    params.native_dir = os.path.abspath(json_object["native_dir"])
    os.chdir(cur_cwd)

    params.sources = json_object["sources"]
    if not type(params.sources) is list:
        raise Exception("Malformed config: sources must be a list")
    for source in params.sources:
        if source not in all_sources:
            raise Exception("Malformed config: invalid source " + source)
    params.ling_types = json_object["ling_types"]
    if not type(params.ling_types) is list:
        raise Exception("Malformed config: ling_types must be a list")
    for ling_type in params.ling_types:
        if ling_type not in all_ling_types:
            raise Exception("Malformed config: invalid ling_type " + ling_type)
    params.msa_types = json_object["msa_types"]
    if not type(params.msa_types) is list:
        raise Exception("Malformed config: msa_types must be a list")
    for msa_type in params.msa_types:
        if msa_type not in all_msa_types:
            raise Exception("Malformed config: invalid msa_type " + msa_type)
    params.partition_types = json_object["partition_types"]
    if not type(params.partition_types) is list:
        raise Exception("Malformed config: partition_types must be a list")
    for partition_type in params.partition_types:
        if partition_type[0] not in ["bin", "multi", "ambig"] or partition_type[1] not in ["MK", "GTR", "BIN"] or partition_type[2] not in [0, 1] or partition_type[3] not in ["2", "x"]:
            raise Exception("Malformed config: invalid partition_type " + partition_type)

    params.glottolog_tree_required = json_object["glottolog_tree_required"]
    if not params.glottolog_tree_required in [0, 1]:
        raise Exception("Malformed config: glottolog_tree_required must be 0 or 1")
    params.flat_paths = json_object["flat_paths"]
    if not params.flat_paths in [0, 1]:
        raise Exception("Malformed config: flat_paths must be 0 or 1")

    try:
        params.download_cutoff = datetime.fromisoformat(json_object["download_cutoff"])
    except:
        raise Exception("Malformed config: download_cutoff must be in isoformat YYYY-MM-DDThh:mm:ss±hh:mm")

def data():
    metadata_path = pb.metadata_path()
    if not os.path.isfile(metadata_path):
        return None
    else:
        return  pd.read_csv(metadata_path, sep = ";", converters = converters)



def write_csv(data_units):
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
        db_df.at[i, "site_group_sizes"] = data.site_group_sizes()
        db_df.at[i, "multistate_ratio"] = data.get_multistate_ratio()
        db_df.at[i, "value_number_counts"] = data.get_value_number_counts()
        db_df.at[i, "value_number_matrix"] = data.get_value_number_matrix()
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
        x = data.max_values_prototype()
        db_df.at[i, "max_values_prototype"] = x
        db_df.at[i, "MULTIx_MK_prototype"] = "MULTI" + str(x) + "_MK"
        db_df.at[i, "MULTIx_GTR_prototype"] = "MULTI" + str(x) + "_GTR"
        db_df.at[i, "COGx"] = "COG" + str(pow(2, x))
        db_df.at[i, "num_discarded_prototype"] = data.num_discarded_prototype()
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

        sampled_msa_paths = []
        for j in range(params.num_samples):
            path = data_unit.sample_path(j)
            if os.path.isfile(path):
                sampled_msa_paths.append(path)
        db_df.at[i, "sampled_msa_paths"] = ""
        db_df.at[i, "sampled_msa_paths"] = sampled_msa_paths

        partition_paths = {}
        for [msa_type, model, gamma, mode] in params.partition_types:
            partition_path = data_unit.partition_path(msa_type, model, gamma, mode)
            if os.path.isfile(partition_path):
                partition_paths[pb.partition_name(msa_type, model, gamma, mode)] = partition_path
        db_df.at[i, "partition_paths"] = ""
        db_df.at[i, "partition_paths"] = partition_paths

        #source information
        db_df.at[i, "primary_source"] = data_unit.primary_source
        db_df.at[i, "native_data_source_sha"] = data_unit.sha
        db_df.at[i, "glottolog_version"] = params.glottolog_version
        db_df.at[i, "download_cutoff"] = params.download_cutoff

    db_df.to_csv(pb.metadata_path(), sep = ";")





def compile():
    pb.clear_data()
    if not (os.path.isdir(pb.domain_path("glottolog")) and os.path.isdir(pb.domain_path("native"))):
        print(colored("Native data does not exist, run download before", "red"))
        return

    if "ambig" in params.msa_types:
        state_encoding.write_charmaps()
    data_units = []
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in os.listdir(os.path.join(pb.domain_path("native"), ds_id)):
            if not source in params.sources:
                continue
            for ling_type in params.ling_types_for_source[source]:
                if not ling_type in params.ling_types:
                    continue
                print(colored("\nChecking data for [ " + ds_id + " " + source + " " + ling_type + " ]", "white"))
                data_units += compile_data_units(ds_id, source, ling_type)
    write_csv(data_units)
    print(colored("Compiling finished. Data saved to " + params.data_dir, "green"))


def compile_data_units(ds_id, source, ling_type):
    handler = native_data.get_handler(ds_id, source)
    data_list = handler.get_data(ling_type)
    data_units = []
    if data_list == []:
        return data_units

    #general
    sha = handler.get_sha()
    primary_source = handler.get_source()

    for data, family in data_list:
        if data is None:
            continue
        data_unit = DataUnit(data)
        data_unit.ds_id = ds_id
        data_unit.source = source
        data_unit.ling_type = ling_type
        data_unit.family = family
        data_unit.sha = sha
        data_unit.primary_source = primary_source


        #tree
        (glottocodes, complete) = handler.get_glottocodes(data.taxon_ids)
        tree_path = data_unit.glottolog_tree_path()
        if complete:
            glottolog_tree = glottolog.get_tree(glottocodes)
            if glottolog_tree is not None:
                pb.mk_file_dir(tree_path)
                glottolog_tree.write(outfile = tree_path, format=9)
        if params.glottolog_tree_required and not os.path.isfile(tree_path):
            print(colored("Missing Glottolog Tree", "red"))
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

        compile_samples(data, data_unit)
        compile_paritions(data, data_unit)
        if "membership" in params.msa_types:
            membership.generate_membership_msa(ds_id, source, ling_type, family)

        data_units.append(data_unit)
        print(colored("Data created for [ " + ds_id + " " + source + " " + ling_type + " " + family + " ]", "green"))

    return data_units


def compile_samples(data, data_unit):
    for i in range(params.num_samples):
        path = data_unit.sample_path(i)
        if i == 0:
            pb.mk_file_dir(path)
        sample = data.get_random_sample(i)
        sample.write_msa(path, "bin")


def compile_paritions(data, data_unit):
    for [msa_type, model, gamma, mode] in params.partition_types:
        partition_path = data_unit.partition_path(msa_type, model, gamma, mode)
        pb.mk_file_dir(partition_path)
        data.write_partitioning(partition_path, msa_type, model, gamma, mode)

def download():
    params.github_token = getpass("Please enter github token: ")
    glottolog.crawl()
    crawler.crawl()
    print(colored("Download finished. Data saved to " + params.native_dir, "green"))



class DataUnit:

    data = None

    ds_id = ""
    source = ""
    ling_type = ""
    family = ""

    primary_source = ""
    sha = ""
    sub_families = []


    def __init__(self, data):
        self.data = data

    def categorical_path(self):
        return pb.categorical_path(self.ds_id, self.source, self.ling_type, self.family)

    def msa_path(self, msa_type):
        return pb.msa_path(self.ds_id, self.source, self.ling_type, self.family, msa_type)

    def glottolog_tree_path(self):
        return pb.glottolog_tree_path(self.ds_id, self.source, self.family)

    def sample_path(self, i):
        return pb.sample_path(self.ds_id, self.source, self.ling_type, self.family, i)

    def partition_path(self, msa_type, model, gamma, mode):
        return pb.partition_path(self.ds_id, self.source, self.ling_type, self.family, msa_type, model, gamma, mode)
