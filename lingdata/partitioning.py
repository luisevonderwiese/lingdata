from lingdata.categorical import CategoricalData
import lingdata.pathbuilder as pb
import os
import pandas as pd


def add_partition_paths(df, model, mode, ambig):
    name = "partition_path_" + pb.partition_name(model, mode, ambig)
    df[name] = ""
    for i, row in df.iterrows():
        data = CategoricalData.from_file(row["categorical_path"])
        if ambig and ("ambig" not in row["msa_paths"] or row["msa_paths"]["ambig"] == ""):
            continue
        partition_path = pb.partition_path(row["ds_id"], row["source"], row["ling_type"], row["family"], model, mode, ambig)
        if os.path.isfile(partition_path):
            df.at[i, name] = partition_path
    return df

def generate_paritions(df, model, mode, ambig):
    for _, row in df.iterrows():
        data = CategoricalData.from_file(row["categorical_path"])
        if ambig and ("ambig" not in row["msa_paths"] or row["msa_paths"]["ambig"] == ""):
            return
        partition_path = pb.partition_path(row["ds_id"], row["source"], row["ling_type"], row["family"], model, mode, ambig)
        pb.mk_file_dir(partition_path)
        data.write_ng_partition(partition_path, model, mode, ambig, pb)
