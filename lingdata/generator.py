import sys
from termcolor import colored

import lingdata.params as params
from lingdata.native_data import CLDFHandler


def set_unlimited_params():
    params.family_split_threshold = sys.maxsize
    params.max_num_taxa = sys.maxsize
    params.max_chars = sys.maxsize


def cldf_to_msa(cldf_path, ling_type, msa_path, msa_type):
    set_unlimited_params()
    handler = CLDFHandler(cldf_path)
    data_list = handler.get_data(ling_type)
    if len(data_list) == 0:
        print(colored("Conversion not possible, check input", "red"))
    assert(len(data_list) == 1)
    data_list[0][0].write_msa(msa_path, msa_type)
    print(colored("Generated " + msa_type + " MSA, saved to " + msa_path, "green"))
