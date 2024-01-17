
from github import Github, UnknownObjectException
import requests
import os
import json
from termcolor import colored

import lingdata.pathbuilder as pb
import lingdata.params as params


cldf_src_file_names = ["cldf/README.md",
                      "cldf/languages.csv",
                      "cldf/values.csv",
                      "cldf/forms.csv",
                      "cldf/cognates.csv"]


cp_user_name = "lingpy"
cp_src_dirs = ["datasets",
                      "trimmed",
                      "data/correspondences"]
cp_dest_file_names = ["dataset.tsv",
                       "trimmed.tsv",
                       "correspondence.tsv"]



def download_file(repo, src_file_name, dest_file_name, sha):
    try:
        url = repo.get_contents(src_file_name, sha).download_url
    except UnknownObjectException:
        return False
    r = requests.get(url, allow_redirects=True)
    open(dest_file_name, 'wb').write(r.content)
    return True




def crawl_cldf():
    github = Github(params.github_token)
    repos = []
    for user in params.source_types["cldf"]:
        if user in params.sources:
            repos += github.get_user(user).get_repos()
    for repo in repos:
        parts = repo.full_name.split("/")
        ds_id = parts[1]
        source = parts[0]
        commits = repo.get_commits(until = params.download_cutoff)
        try:
            sha = commits[0].sha
        except: #No older commits, repo did not exist at cutoff date
            print(colored(ds_id + " from " + source +  " skipped", "yellow"))
            continue
        download_dir = pb.source_path("native", ds_id, source)
        dest_file_names = [os.path.join(download_dir, file_name.split("/")[-1]) for file_name in cldf_src_file_names]
        if os.path.isdir(download_dir):
            meta_path = os.path.join(download_dir, "meta.json")
            if os.path.isfile(meta_path):
                with open(meta_path, 'r') as openfile:
                    json_data = json.load(openfile)
                if json_data["sha"] == sha:
                    print(colored(ds_id + " from " + source +  " up to date", "yellow"))
                    continue
        pb.mk_this_dir(download_dir)
        meta_dict = {"sha" : sha}
        with open(os.path.join(download_dir, "meta.json"), 'w+') as outfile:
            json.dump(meta_dict, outfile)
        try:
            for (i, src_file_name) in enumerate(cldf_src_file_names):
                download_file(repo, src_file_name, dest_file_names[i], sha)
            print(colored(ds_id + " from " + source +  " downloaded", "green"))
        except:
            print(colored(ds_id + " from " + source +  ": error occured", "red"))
            pb.rm_this_dir(download_dir)



def crawl_cp():
    source = params.source_types["correspondence"][0]
    if not source in params.sources:
        return
    repo = Github(params.github_token).get_user(cp_user_name).get_repo(source)
    commits = repo.get_commits(until = params.download_cutoff)
    try:
        sha = commits[0].sha
    except: #No older commits, repo did not exist at cutoff date
        print(colored("correspondence-pattern-data skipped", "yellow"))
        return
    #all files from same repo, so check only at first file
    some_ds_id = repo.get_contents(cp_src_dirs[0], sha)[0].path.split("/")[-1].split(".")[0]
    meta_path = os.path.join(pb.source_path("native", some_ds_id, source), "meta.json")
    if os.path.isfile(meta_path):
        with open(meta_path, 'r') as openfile:
            json_data = json.load(openfile)
        if json_data["sha"] == sha:
            print(colored("correspondence-pattern-data up to date", "yellow"))
            return
    meta_dict = {"sha" : sha}

    for (i, src_dir) in enumerate(cp_src_dirs):
        repo_contents = repo.get_contents(src_dir, sha)
        for content_file in repo_contents:
            ds_id = content_file.path.split("/")[-1].split(".")[0]
            source_path = pb.source_path("native", ds_id, source)
            dest_file_name = os.path.join(source_path, cp_dest_file_names[i])
            pb.mk_file_dir(dest_file_name)
            try:
                download_file(repo, content_file.path, dest_file_name, sha)
                print(colored(ds_id + " from correspondence-pattern-data downloaded", "green"))
            except:
                print(colored(ds_id + " from correspondence-pattern-data: error occured", "red"))
                pb.rm_this_dir(source_path)
                return
            if i == 0:
                meta_path = os.path.join(source_path, "meta.json")
                with open(meta_path, 'w+') as outfile:
                    json.dump(meta_dict, outfile)

def crawl():
    crawl_cldf()
    crawl_cp()
