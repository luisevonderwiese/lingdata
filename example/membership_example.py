import lingdata.database as database

config_path = "membership_lingdata_config.json"
database.read_config(config_path)
#database.download()
database.compile()
df = database.data()
print(df)
