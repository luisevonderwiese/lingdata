import sys
sys.path.append("..")
import lingdata.database as database
import lingdata.params as params
import lingdata.sampling as sampling
import lingdata.partitioning as partitioning
import os
import pandas as pd
import json

pd.set_option('display.max_rows', None)

config_path = "lingdata_example_config.json"
database.read_config(config_path)
database.update_native()
database.generate_data()
df = database.data()


with open(config_path, 'r') as openfile:
    json_object = json.load(openfile)
sampling.generate_samples(df, 10)
df = sampling.add_sampled_msa_paths(df, 10)


partition_types =  [
                      ["MK", "x", 0],
                      ["MK", "2", 0],
                      ["GTR", "x", 0],
                      ["GTR", "2", 0],
                      #["MK", "x", 1],
                      #["MK", "2", 1],
                      #["GTR", "x", 1],
                      #["GTR", "2", 1]
                      ]

for [model, mode, ambig] in partition_types:
    partitioning.generate_paritions(df, model, mode, ambig)
    df = partitioning.add_partition_paths(df, model, mode, ambig)
