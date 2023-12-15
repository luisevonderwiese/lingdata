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
#database.update_native()
#database.generate_data()
df = database.data()
print(df)
