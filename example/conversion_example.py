import lingdata.generator as generator

cldf_path = "../conversion_example_data/cldf"
msa_path = "../conversion_example_data/msa/bin.phy"
ling_type = "cognate"
msa_type = "bin"

generator.cldf_to_msa(cldf_path, ling_type, msa_path, msa_type)
