import os
import shutil
import lingdata.params as params

native_domains = ["native", "glottolog"]
def mk_file_dir(path):
    super_dir = "/".join(path.split("/")[:-1])
    mk_this_dir(super_dir)

def mk_this_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def rm_this_dir(path):
    shutil.rmtree(path)


def domain_path(domain):
    if domain in native_domains:
        return os.path.join(params.native_dir, domain)
    else:
        return os.path.join(params.data_dir, domain)


def dataset_path(domain, ds_id):
    return os.path.join(domain_path(domain), ds_id)

def source_path(domain, ds_id, source):
    if not domain in native_domains and params.flat_paths:
        return "_".join([dataset_path(domain, ds_id), source])
    else:
        return os.path.join(dataset_path(domain, ds_id), source)

def type_path(domain, ds_id, source, ling_type):
    if params.flat_paths:
        return "_".join([source_path(domain, ds_id, source), ling_type])
    else:
        return os.path.join(source_path(domain, ds_id, source), ling_type)

def family_path(domain, ds_id, source, ling_type, family):
    if params.flat_paths:
        return "_".join([type_path(domain, ds_id, source, ling_type), family])
    else:
        return os.path.join(type_path(domain, ds_id, source, ling_type), family)

def metadata_path():
    return os.path.join(params.data_dir, "lingdata.csv")

def glottolog_tree_path(ds_id, source, family):
    if params.flat_paths:
        return os.path.join("_".join([source_path("glottolog_trees", ds_id, source), family]), "glottolog.tre")
    else:
        return os.path.join(source_path("glottolog_trees", ds_id, source), family, "glottolog.tre")

def categorical_path(ds_id, source, ling_type, family):
    return os.path.join(family_path("categorical", ds_id, source, ling_type, family), "categorical.csv")

def msa_path(ds_id, source, ling_type, family, msa_type):
    if msa_type == "catg_bin":
        name = "bin.catg"
    elif msa_type == "catg_multi":
        name = "multi.catg"
    else:
        name = msa_type + ".phy"
    return os.path.join(family_path("msa", ds_id, source, ling_type, family), name)

def partition_path(ds_id, source, ling_type, family, multi_model, mode, ambig):
    return os.path.join(family_path("partition", ds_id, source, ling_type, family), partition_name(multi_model, mode, ambig) + ".part")

def partition_name(multi_model, mode, ambig):
    name = multi_model + mode
    if ambig:
        name += "+M"
    return name

def sample_path(ds_id, source, ling_type, family, idx):
    return os.path.join(family_path("sampled", ds_id, source, ling_type, family), "sampled" + str(idx) + "_bin.phy")


def charmap_path(x):
    return os.path.join(domain_path("charmaps"), "charmap_" + str(x) + ".txt")


def clear_data():
    if os.path.isdir(params.data_dir):
        shutil.rmtree(params.data_dir)
