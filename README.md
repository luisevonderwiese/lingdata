# lingdata
This is a python package for handling linguistic data for phylogenetic inference

The package allows it to download linguistic data from different sources automatically and to compile a database, containing different types of multiple sequence alignments (MSAs) / character matrices and reference trees extracted from glottolog.

The package can also be used to generate different types of MSAs / character matrices for linguistic data provided locally in CLDF Format.


## Installation
```
pip install .
```
## Usage within Python
See `example/lingdata_example.py`

## Usage from Command Line

### Database
```
lingdata --download -c example/lingdata_example_config.json
lingdata --compile -c example/lingdata_example_config.json  

```


### Generating MSAs / character matrices
```
lingdata --generate -i conversion_example_data/cldf/ -l cognate -o conversion_example_data/msa/bin.phy  -m bin
```




## Config
For an example see `example/lingdata_example_config.json`
| parameter name | explanation |
| --- | --- |
| `max_num_taxa` | maximum number of languages such that the dataset is added to the compiled database |
| `max_num_chars` | maximum number of concepts such that the dataset is added to the compiled database |
| `family_split_threshold` | maximum number of languages such that the dataset is added to the compiled database without being split in subfamilies |
| `num_samples` | number of sampled binary MSAs in the compiled database (see Paper) |
| `data_dir` | directory where generated data is stored <br>(absolute or relative to location of config file)|
| `native_dir` | directory where data downloaded from sources is stored <br>(absolute or relative to location of config file)|
| `sources` | sources from which data is added to the compiled database <br>Supported sources: `"lexibank"` (<https://github.com/lexibank>), `"SequenceComparison"` (<https://github.com/SequenceComparison>), `"correspondence-pattern-data"` (<https://github.com/lingpy/correspondence-pattern-data>)<sup>1</sup> |
| `ling_types` | Datasets of the provided linguistic data type are added to the compiled database <br>Supported types: `"cognate"`, `"structural"` (morpho-syntactic or morpho-phonological data), `"correspondence"` (sound correspondence patterns) |
| `msa_types` | MSA (character matrix) types contained in the compiled database <br>Supported types: `"bin"` (binary), `"multi"` (multi-valued), `"catg_bin"` (probabilistic binary), `"catg_multi"` (probabilistic multi-valued), `"ambig"` (with user defined state encoding)<sup>2</sup> |
| `partition_types` | [Partitionings]() contained in the compiled database <br>Partition type specified as [msa_type, multi_model, gamma, mode] where <br>msa_type: valid msa_type (see above) <br>multi_model: `"MK"`, `"GTR"` ([details](https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model)) <br>gamma: `0`, `1` ([details](https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model)) <br>mode: `"2"`, `"x"` (`"2"`: two partitions only (characteristics with <= 2 / > 2 values), "x"`: One partition for each number of values per characteristic)|
| `glottolog_tree_required` | `1`: only datasets with existing glottolog tree must be added to the compiled database <br>`0`: datasets are added to the compiled database no matter whether a glottolog tree exists or not |
| flat_paths | `1`: data is stored in a directory structure with one level only <br>`0`: data is stored in hierarchical structure with one level for each ds_id, source, ling_type, family |

<sup>1</sup>For downloading the native data from Github, a Github Token is required. The repo correspondence-pattern-data is private, so this source can only be used if the Github user has access to it.

<sup>2</sup>Not every type can be created for every data set. Details [here](https://github.com/amkozlov/raxml-ng/wiki/Input-data#multiple-sequence-alignment)


## File Structure of the Compiled Database
Each of the listed directories (except from `native_dir/glottolog/` and `data_dir/charmaps/`) contains a subdirectory with the respective files for each dataset. All paths are provided in the respective columns of `data_dir/lingdata.csv`.
| directory | contains |
| --- | --- |
| `native_dir/glottolog/` | List of languages and full gold standard tree from glottolog |
| `native_dir/native/` | Data downloaded from sources (in CLDF format mainly) |
| `data_dir/lingdata.csv` | Actual database file containing metadata and paths |
| `data_dir/categorical/` | Categorical matrices (required as an intermediate step for creating MSAs / character matrices)
| `data_dir/glottolog_trees/` | Glottolog reference trees pruned from the full tree |
| `data_dir/msa/` | Different types of MSAs / character matrices |
| `data_dir/sampled/` | Sampled binary MSAs |
| `data_dir/partitionings/` | Paritition files |
| `data_dir/charmaps/` | Charmaps for custom state encoding |
