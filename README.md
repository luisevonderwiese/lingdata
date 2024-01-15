# lingdata
This is a python package for handling linguistic data for phylogenetic inference
## Installation
```
pip install .
```
## Usage
See `example/lingdata_example.py`

## Config
| parameter name | explanation |
|--- | --- |
| max_num_taxa | maximum number of languages such that the dataset is added to the compiled database |
| max_num_chars | maximum number of concepts such that the dataset is added to the compiled database |
| family_split_threshold | maximum number of languages such that the dataset is added to the compiled database without being split in subfamilies |
| num_samples | number of sampled binary MSAs (see Paper) |
| data_dir | directory where generated data is stored |
| native_dir | directory where data downloaded from sources is stored |
| sources | supported sources: "lexibank", "SequenceComparison", "correspondence-pattern-data" |
| ling_types | supported types: "cognate", "structural", "correspondence" (sound correspondence patterns) |
| msa_types | supported types: "bin", "multi", "catg_bin", "catg_multi", "ambig" |
| partition_types |  a list of partition types specified as [msa_type, multi_model, gamma, mode] where <br>msa_type: valid msa_type (see above) <br>multi_model: in ["MK", "GTR"] <br>gamma in [0, 1] <br>mode in ["2", "x"] |
| glottolog_tree_required | 1: only datasets with existing glottolog tree must be added to the compiled database <br>0: datasets are added to the compiled database no matter whether a glottolog tree exists or not |
| flat_paths | 1: data is stored in a directory structure with one level only <br>0: data is stored in hierarchical structure with one level for each ds_id, source, ling_type, family |
