max_num_taxa = float("nan")
max_num_chars = float("nan")
family_split_threshold = float("nan")
num_samples = float("nan")

data_dir = ""
native_dir = ""

sources = []
ling_types = []

msa_types = []
partition_types = []


glottolog_tree_required = False

flat_paths = False


github_token = ""

glottolog_version = "v4.8"
glottolog_tree_url = "https://cdstar.eva.mpg.de//bitstreams/EAEA0-B701-6328-C3E3-0/tree_glottolog_newick.txt"

max_max_values = 64 #restriction of RAxML-NG
min_num_taxa = 4 #restriction of RAxML-NG

source_types = {"cldf" : ["lexibank", "SequenceComparison"], "correspondence" : ["correspondence-pattern-data"]}

ling_types_for_source = {"lexibank" : ["cognate", "structural"], "SequenceComparison": ["cognate", "structural"], "correspondence-pattern-data" : ["cognate", "correspondence"]}
