import os
import requests
import copy
import re
import pandas as pd
from ete3 import Tree
import numpy as np
from github import Github, UnknownObjectException
from termcolor import colored
import json

import lingdata.params as params
import lingdata.pathbuilder as pb


def raw_tree_path():
    return os.path.join(pb.domain_path("glottolog"), "tree_glottolog_newick.txt")
def full_tree_path():
    return os.path.join(pb.domain_path("glottolog"), "glottolog.tre")


cldf_repo_name = "glottolog/glottolog-cldf"
required_files = ["languages.csv"]

family_dict = None
full_tree = None


def download_file(repo, file_name, sha):
    download_path = os.path.join(pb.domain_path("glottolog"), file_name)
    if os.path.isfile(download_path):
        return True
    try:
        url = repo.get_contents(os.path.join("cldf", file_name), sha).download_url
    except UnknownObjectException:
        return False
    r = requests.get(url, allow_redirects=True)
    open(download_path, 'wb').write(r.content)
    return True


def crawl():
    glottolog_path = pb.domain_path("glottolog")
    github = Github(params.github_token)
    cldf_repo = github.get_repo(cldf_repo_name)
    tag = next((x for x in cldf_repo.get_tags() if x.name == params.glottolog_version))
    sha = tag.commit.sha
    if os.path.exists(glottolog_path):
        meta_path = os.path.join(glottolog_path, "meta.json")
        if os.path.isfile(meta_path):
            with open(meta_path, 'r') as openfile:
                json_data = json.load(openfile)
            if json_data["sha"] == sha:
                print(colored("Glottolog files up to date", "yellow"))
            return
            pb.rm_this_dir(glottolog_path)
    pb.mk_this_dir(glottolog_path)
    meta_dict = {"sha" : sha}
    with open(os.path.join(glottolog_path, "meta.json"), 'w+') as outfile:
        json.dump(meta_dict, outfile)

    for file_name in required_files:
        try:
            download_file(cldf_repo, file_name, sha)
            print(colored("Glottolog " + file_name + " downloaded", "green"))
        except:
            print(colored("Glottolog " + file_name + ": error occured", "red"))
            pb.rm_this_dir(glottolog_path)
            return
    try:
        r = requests.get(params.glottolog_tree_url, allow_redirects=True)
        open(raw_tree_path(), 'wb').write(r.content)
        print(colored("Glottolog tree downloaded", "green"))
    except:
        print(colored("Glottolog tree: error occured", "red"))
        pb.rm_this_dir(glottolog_path)
        return
    try:
        extract_full_glottolog_tree()
        print(colored("Glottolog tree extracted", "green"))
    except:
        print(colored("Glottolog tree extraction: error occured", "red"))
        pb.rm_this_dir(glottolog_path)
        return

def extract_full_glottolog_tree():
    #code adapted from gerhard jaeger
    with open(raw_tree_path()) as f:
        raw = f.readlines()
    trees = []
    # each line is a tree. bring in proper format and read with ete3
    for i, ln in enumerate(raw):
        ln = ln.strip()
        ln = re.sub(r"\'[A-Z][^[]*\[", "[", ln)
        ln = re.sub(r"\][^']*\'", "]", ln)
        ln = re.sub(r"\[|\]", "", ln)
        ln = ln.replace(":1", "")
        trees.append(Tree(ln, format=1))
    # place all trees below a single root
    glot = Tree()
    for t in trees:
        glot.add_child(t)

    #insert missing, i.e. isolated languages below the root
    tTaxa = [nd.name for nd in glot.traverse() if nd.name != '']
    fn = os.path.join(pb.domain_path("glottolog"), "languages.csv")
    languages = pd.read_csv(fn)
    gTaxa = languages.Glottocode.values
    for taxon in gTaxa:
        if taxon not in tTaxa:
            glot.add_child(name=taxon)

    #if there is a inner node with a name (i.e. corresponds to a language),the name of this node is removed
    # and a child(i.e.leaf) with this name is inserted
    nonLeaves = [nd.name for nd in glot.traverse() if nd.name != '' and not nd.is_leaf()]
    for i, nm in enumerate(nonLeaves):
        nd = glot & nm
        nd.name = ''
        nd.add_child(name=nm)

    # only keep languages which are listed in languages.csv
    gTaxa = np.intersect1d(gTaxa, glot.get_leaf_names())
    glot.prune([glot&x for x in gTaxa])

    glot.write(outfile = full_tree_path(), format=9)
    global full_tree
    full_tree = glot


def load_families():
    global family_dict
    if family_dict is not None:
        return True
    languages_path = os.path.join(pb.domain_path("glottolog"), "languages.csv")
    if not os.path.isfile(languages_path):
        return False
    languages_df = pd.read_csv(languages_path)
    family_dict = {}
    for index, row in languages_df.iterrows():
        glottocode = row["Glottocode"]
        family_id = str(row["Family_ID"])
        if family_id == "nan":
            family_dict[glottocode] = "ISOLATE"
        else:
            family_dict[glottocode] = family_id
    return True

def load_full_tree():
    global full_tree
    if full_tree is not None:
        return True
    if not os.path.isfile(full_tree_path()):
        return False
    full_tree = Tree(full_tree_path(), format=9)
    return True



def get_families(glottocodes):
    r = load_families()
    family_ids = set()
    if not r:
        return family_ids
    for glottocode in glottocodes:
        if glottocode not in family_dict:
            family_ids.add("UNKNOWN")
            continue
        family_ids.add(family_dict[glottocode])
    return family_ids

def split_families(glottocodes):
    r = load_families()
    families = {}
    if not r:
        return families
    for glottocode, lang_ids in glottocodes.items():
        if glottocode not in family_dict:
            family_id = "UNKNOWN"
        else:
            family_id = family_dict[glottocode]
        if family_id in families:
            families[family_id] += lang_ids
        else:
            families[family_id] = lang_ids
    return families



def get_tree(glottocodes):
    r = load_full_tree()
    if not r:
        return None
    tree = copy.deepcopy(full_tree)
    try:
        tree.prune([tree&glottocode for glottocode in glottocodes])
    except: #node not found due to wrong / deprecated glottocodes
        return None
    for leaf in tree.iter_leaves():
        leaf.add_features(new = False)
    for leaf in tree.iter_leaves():
        if leaf.new:
            continue
        ids = glottocodes[leaf.name]
        if len(ids) == 1:
            leaf.name = ids[0]
        else:
            leaf.name = ""
            for i in ids:
                leaf.add_child(name = i)
            for child in leaf.children:
                child.add_features(new = True)
    return tree
