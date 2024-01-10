# lingdata
This is a python package for handling linguistic data for phylogenetic inference
## Installation
```
pip install .
```
## Usage
See `example/lingdata_example.py`

## Config

| max_num_taxa | maximum number of languages such that the dataset is added to the compiled database |
| max_num_chars | maximum number of concepts such that the dataset is added to the compiled database |
| family_split_threshold | maximum number of languages such that the dataset is added to the compiled database without being split in subfamilies |
| num_samples | number of sampled binary MSAs (see Paper) |
| data_dir | directory where generated data is stored |
| native_dir | directory where data downloaded from sources is stored |
| sources | lexibank", "SequenceComparison", "correspondence-pattern-data" |
| ling_types | cognate", "structural", "correspondence" |
| msa_types | bin", "multi", "catg_bin", "catg_multi", "ambig" |
| partition_types |  a list of partition types specified in the structure [msa_type, multi_model, gamma, mode] where msa_type must be a valid msa_type (see above), multi_model in ["MK", "GTR"], gamma in [0, 1], mode in ["2", "x"] |
| glottolog_tree_required | 1 if only datasets with existing glottolog tree must be added to the compiled database, 0 otherwise |
| flat_paths | 1 - data stored in a directory structure with one level only, 0 - data stored in hierarchical structure with one level for each ds_id, source, ling_type, family |
