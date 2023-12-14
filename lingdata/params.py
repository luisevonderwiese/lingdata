size_limit = float("nan")
family_split_threshold = float("nan")

data_dir = ""
native_dir = ""

sources = []
ling_types = []

msa_types = []
partition_types = []


glottolog_tree_required = False

flat_paths = False


token = "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" #for github

glottolog_version = "v4.8"
glottolog_tree_url = "https://cdstar.eva.mpg.de//bitstreams/EAEA0-B701-6328-C3E3-0/tree_glottolog_newick.txt"

max_max_values = 64 #do not use datasets in which characteristics have too many values...will be problematic

source_types = {"cldf" : ["lexibank", "SequenceComparison"], "correspondence" : ["correspondence-pattern-data"]}
