import lingdata.database as database

config_path = "lingdata_example_config.json"
database.read_config(config_path)
#database.update_native()
database.generate_data()
df = database.data()
print(df)
