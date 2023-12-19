import sys
sys.path.append("..")
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
import os
import shutil

from lingdata.native_data import CLDFHandler
from lingdata.categorical import CategoricalData
import lingdata.glottolog as glottolog
import lingdata.database as database
import lingdata.params as params
import lingdata.pathbuilder as pb
import lingdata.native_data as native_data


from oldcldf import OldCLDFHandler
import correspondences_jaeger
import partitioning_reference





def equal_aligns(align1, align2):
    #align1 = replace_indef_zero(align1)
    num_sites1 = align1.get_alignment_length()
    num_sites2 = align2.get_alignment_length()
    if num_sites1 != num_sites2:
        print(num_sites1)
        print(num_sites2)
        print("!!! Different number of sites")
        return False
    #align1 = taxa_name_unification(align1)
    #align2 = taxa_name_unification(align2)
    taxa1 = set([align1[i].id for i in range(len(align1))])
    taxa2 = set([align2[i].id for i in range(len(align2))])
    common_taxa = taxa1.intersection(taxa2)
    taxa1_only = taxa1.difference(common_taxa)
    taxa2_only = taxa2.difference(common_taxa)
    if (len(taxa1_only) != 0 or len(taxa2_only) != 0):
        print(taxa1_only)
        print(taxa2_only)
        print("!!! Different Taxa")
        return False
    align1.sort(key=lambda x: x.id)
    align2.sort(key=lambda x: x.id)
    columns1 = [align1[:, i] for i in range(num_sites1)]
    columns2 = [align2[:, i] for i in range(num_sites2)]

    L1 = [ (columns1[i],i) for i in range(len(columns1)) ]
    L1.sort()
    columns1, permutation1 = zip(*L1)

    L2 = [ (columns2[i],i) for i in range(len(columns2)) ]
    L2.sort()
    columns2, permutation2 = zip(*L2)

    j = 0
    i = 0
    while i < len(columns1) and j < len(columns2):
        if columns1[i] == columns2[j]:
            i+=1
            j+=1
            continue
        if columns1[i] < columns2[j]:
            print("Unmatch 1")
            print(columns1[i])
            #print(permutation1[i])
            i+=1
            return False
        if columns1[i] > columns2[j]:
            print("Unmatch 2")
            print(columns2[j])
            #print(permutation2[j])
            j+=1
            return False
    return True

def equal_files(path1, path2):
    f1 = open(path1, "r")
    f2 = open(path2, "r")
    f1_data = f1.readlines()
    f2_data = f2.readlines()
    if len(f1_data) != len(f2_data):
        print("Different number of lines")
        f1.close()
        f2.close()
        return False
    for (i, l1) in enumerate(f1_data):
        if l1 != f2_data[i]:
            print("Difference in line " + str(i))
            print(l1)
            print(f2_data[i])
            f1.close()
            f2.close()
            return False
    f1.close()
    f2.close()
    return True

def is_sample(this, other):
    if not other.is_single_state():
        print("Not single state")
        return False
    if this.char_ids != other.char_ids:
        print("char_ids")
        print(this.char_ids)
        print(other.char_ids)
        return False
    if this.taxon_ids != other.taxon_ids:
        print("taxon_ids")
        print(this.taxon_ids)
        print(other.taxon_ids)
        return False
    for char_idx in range(this.num_chars()):
        for taxon_idx in range(this.num_taxa()):
            if len(other.matrix[char_idx][taxon_idx]) == 0:
                if len(this.matrix[char_idx][taxon_idx]) != 0:
                    return False
                continue
            sampled_value = other.matrix[char_idx][taxon_idx][0]
            if sampled_value not in this.matrix[char_idx][taxon_idx]:
                print(sampled_value + " not in " + str(this.matrix[char_idx][taxon_idx]))
                return False
    return True


def test_glottolog():
    glottocodes = {"siaw1243" : ["siaw1243"],
                   "amto1250" : ["amto1250"],
                   "khak1248" : ["khak1248"],
                   "west2402" : ["west2402"],
                   "ainu1251" : ["ainu1251"],
                   "nhen1239" : ["nhen1239"],
                   "tupi1273" : ["tupi1273"],
                   "hier1240" : ["hier1240"]}

    glottolog.get_families(glottocodes)
    glottolog.split_families(glottocodes)
    glottolog.get_tree(glottocodes)
    print("========== GLOTTOLOG FUNCTIONAL TEST PASSED")


def test_categorical_read_write():

    test_dir = "test_temp/"
    path = os.path.join(test_dir, "test.csv")
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.sources:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    data.write(path)
                    data2 = CategoricalData.from_file(path)
                    assert(data == data2)
    shutil.rmtree(test_dir)
    print("========== CATEGROICAL READ WRITE CORRECTNESS TEST PASSED")



def test_df_creation():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    if family != "full":
                        continue
                    source_path = pb.source_path("native", ds_id, source)
                    old_cldf = OldCLDFHandler(source_path)
                    old_data = CategoricalData.from_list_df(old_cldf.get_df(ling_type))
                    assert(old_data == data)
    print("========== CLDF CATEGORICAL CORRECTNESS TEST PASSED")


def test_df_creation_correspondence():
    test_dir = "test_temp/"
    path = os.path.join(test_dir, "test.csv")
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["correspondence"]:
            ling_type = "correspondence"
            handler = native_data.get_handler(ds_id, source)
            data_list = handler.get_data(ling_type)
            for data, family in data_list:
                if data is None:
                    continue
                if family != "full":
                    continue
                align1 = data.get_msa("bin")
                align2 = correspondences_jaeger.bin_align(os.path.join(pb.source_path("native", ds_id, source), "correspondence.tsv"))
                assert(equal_aligns(align1, align2))
    print("========== CORRESPONDENCES CATEGORICAL CORRECTNESS TEST PASSED")


def test_bin_align():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    if family != "full":
                        continue
                    align1 = data.get_msa("bin")
                    old_cldf = OldCLDFHandler(pb.source_path("native", ds_id, source))
                    align2 = old_cldf.get_bin_align(ling_type)
                    assert(equal_aligns(align1, align2))
    print("========== BINARY MSA CORRECTNESS TEST PASSED")


def test_multi_align():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    if family != "full":
                        continue
                    align1 = data.get_msa("multi")
                    old_cldf = OldCLDFHandler(pb.source_path("native", ds_id, source))
                    r, align2 = old_cldf.get_multi_align(ling_type)
                    if align1 is None:
                        assert(r["multi_cells"] > 0 or r["high_state_chars"] > 0 or r["max_states"] < 3)
                    else:
                        assert(align2 is not None)
                        assert(r["max_states"] > 2)
                        assert(r["multi_cells"] == 0)
                        assert(r["high_state_chars"] == 0)
                        assert(equal_aligns(align1, align2))
    print("========== MULTI-VALUE MSA CORRECTNESS TEST PASSED")


def test_ambig_align():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    align = data.get_msa("ambig")
                    #only a functional test, no correctness
    print("========== AMBIG MSA FUNCTIONAL TEST PASSED")



def test_catg():
    test_dir = "test_temp/"
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    old_path = os.path.join(test_dir, "old.csv")
    new_path = os.path.join(test_dir, "new.csv")
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    if family != "full":
                        continue
                    data.write_catg_msa(new_path)
                    old_cldf = OldCLDFHandler(pb.source_path("native", ds_id, source))
                    old_cldf.write_catg(ling_type, old_path)
                    assert(equal_files(new_path, old_path))
    shutil.rmtree(test_dir)
    print("========== BINARY CATG MSA CORRECTNESS TEST PASSED")


def test_multi_catg():
    test_dir = "test_temp/"
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    old_path = os.path.join(test_dir, "old.csv")
    new_path = os.path.join(test_dir, "new.csv")
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.source_types["cldf"]:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    if family != "full":
                        continue
                    data.write_multi_catg_msa(os.path.join(test_dir, ds_id + ".csv")) #only a functional test, no correctness
    shutil.rmtree(test_dir)
    print("========== MULTI-VALUE CATG MSA FUNCTIONAL TEST PASSED")


def test_sample():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.sources:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    for i in range(params.num_samples):
                        sample = data.get_random_sample(i)
                        assert(is_sample(data, sample))
    print("========== SAMPLING CORRECTNESS TEST PASSED")

def test_family_splitting():
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.sources:
            for ling_type in params.ling_types:
                params.family_split_threshold = 5
                handler = native_data.get_handler(ds_id, source)
                split_data_list = handler.get_data(ling_type)
                params.family_split_threshold = 10000000000
                unsplit_data_list = handler.get_data(ling_type)
                if len(unsplit_data_list) == 0:
                    continue
                if len(split_data_list) == 0 or (len(split_data_list) == 1 and split_data_list[0][1] == "full"):
                    continue
                assert(len(unsplit_data_list) == 1 and unsplit_data_list[0][1] == "full")
                full_data = unsplit_data_list[0][0]
                full_subfamilies = glottolog.get_families(handler.get_glottocodes(full_data.taxon_ids)[0])
                for data, family in split_data_list:
                    if data is None:
                        continue
                    assert(family in full_subfamilies)
                    glottocodes, complete = handler.get_glottocodes(data.taxon_ids)
                    assert(complete)
                    for glottocode in glottocodes.keys():
                        assert(glottolog.family_dict[glottocode] == family)
                    for (full_char_idx, char_id) in enumerate(full_data.char_ids):
                        if char_id not in data.char_ids:
                            for taxon_id in data.taxon_ids:
                                full_taxon_idx = full_data.taxon_ids.index(taxon_id)
                                assert(full_data.matrix[full_char_idx][full_taxon_idx] == [])
                        else:
                            char_idx = data.char_ids.index(char_id)
                            for (taxon_idx, taxon_id) in enumerate(data.taxon_ids):
                                full_taxon_idx = full_data.taxon_ids.index(taxon_id)
                                assert(data.matrix[char_idx][taxon_idx] == full_data.matrix[full_char_idx][full_taxon_idx])
    print("========== FAMILIY SPLITTING CORRECTNESS TEST PASSED")

def test_partitioning():
    test_dir = "test_temp/"
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    old_path = os.path.join(test_dir, "old.csv")
    msa_path = os.path.join(test_dir, "msa.phy")
    new_path = os.path.join(test_dir, "new.csv")
    for ds_id in os.listdir(pb.domain_path("native")):
        for source in params.sources:
            for ling_type in params.ling_types:
                handler = native_data.get_handler(ds_id, source)
                data_list = handler.get_data(ling_type)
                for data, family in data_list:
                    if data is None:
                        continue
                    res = data.write_msa(msa_path, "multi")
                    if res:
                        data.write_partitioning(new_path, "multi", "GTR", False, "2")
                        partitioning_reference.create_ng_partition(msa_path, old_path, multi_model = "GTR", complete = False)
                        if os.path.isfile(old_path):
                            assert(equal_files(new_path, old_path))
                            os.remove(old_path)
                        data.write_partitioning(new_path, "multi", "GTR", False, "x")
                        partitioning_reference.create_ng_partition(msa_path, old_path, multi_model = "GTR", complete = True)
                        if os.path.isfile(old_path):
                            assert(equal_files(new_path, old_path))
                            os.remove(old_path)
    shutil.rmtree(test_dir)
    print("========== PARTITIONING (CORRECTNESS) TEST PASSED")





config_path = "lingdata_test_config.json"
database.read_config(config_path)



test_glottolog()
test_categorical_read_write()
test_df_creation()
test_df_creation_correspondence()
test_bin_align()
test_multi_align()
test_ambig_align()
test_catg()
test_multi_catg()
test_sample()
test_family_splitting()
test_partitioning()
