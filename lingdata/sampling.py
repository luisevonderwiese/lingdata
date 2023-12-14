from lingdata.categorical import CategoricalData
import lingdata.pathbuilder as pb
import os
import pandas as pd


def add_sampled_msa_paths(df, num_samples):
    df["sampled_msa_paths"] = ""
    for k, row in df.iterrows():
        df.at[k, "sampled_msa_paths"] = []
        data = CategoricalData.from_file(row["categorical_path"])
        for i in range(num_samples):
            path = pb.sample_path(row["ds_id"], row["source"], row["ling_type"], row["family"], i)
            if os.path.isfile(path):
                df.at[k, "sampled_msa_paths"].append(path)
    return df

def generate_samples(df, num_samples):
    for _, row in df.iterrows():
        data = CategoricalData.from_file(row["categorical_path"])
        for i in range(num_samples):
            path = pb.sample_path(row["ds_id"], row["source"], row["ling_type"], row["family"], i)
            if i == 0:
                pb.mk_file_dir(path)
            sampled = data.get_random_sample(i)
            sampled.write_msa(path, "bin")
