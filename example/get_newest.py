from datetime import datetime
import json
import os
import pytz

import lingdata.pathbuilder as pb
import lingdata.database as database

config_path = "../example/lingdata_example_config.json"

database.read_config(config_path)
newest_date = pytz.utc.localize(datetime.min)
for ds_id in os.listdir(pb.domain_path("native")):
    for source in os.listdir(os.path.join(pb.domain_path("native"), ds_id)):
        meta_path = os.path.join(pb.domain_path("native"), ds_id, source, "meta.json")
        with open(meta_path, 'r') as openfile:
            json_data = json.load(openfile)
        newest_date = max(newest_date, datetime.fromisoformat(json_data["updated_at"]))
print(newest_date.isoformat())
