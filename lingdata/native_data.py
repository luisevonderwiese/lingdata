import os
import json
import pandas as pd
from termcolor import colored

from lingdata.categorical import CategoricalData
import lingdata.params as params
import lingdata.glottolog as glottolog
import lingdata.pathbuilder as pb


def drop_columns_except(df, relevant_columns):
    new_df = df
    for column in df.columns:
        if not column in relevant_columns:
            new_df = new_df.drop(column, axis=1)
    return new_df

class NativeData:

    def get_data(self, families):
        if self.df is None:
            print(colored("Native data incomplete", "red"))
            return []
        c = self.num_chars()
        t = self.num_taxa()
        if c > params.max_num_chars:
            print(colored("Too many chars: " +  str(c), "red"))
            return []
        if families != {} and t > params.family_split_threshold:
            print(colored("Splitting families (" +  str(t) +  ")", "white"))
            return self.split_data(families)
        else:
            if t > params.max_num_taxa:
                print(colored("Too many taxa: " +  str(t), "red"))
                return []
            if t < params.min_num_taxa:
                print(colored("Too few taxa: " +  str(t), "red"))
                return []
            return self.full_data()

    def num_taxa(self):
        return len(self.df['Language_ID'].unique())

    def num_chars(self):
        return len(self.df['Char_ID'].unique())

    def split_data(self, families):
        cds = []
        for family_id, lang_ids in families.items():
            sub_df =  self.df[self.df['Language_ID'].isin(lang_ids)]
            num_taxa = len(sub_df['Language_ID'].unique())
            if num_taxa > params.max_num_taxa:
                print(colored(family_id + " - Too many taxa: " +  str(num_taxa), "red"))
                continue
            if num_taxa < params.min_num_taxa:
                print(colored(family_id + " - Too few taxa: " +  str(num_taxa), "red"))
                continue
            cds.append((CategoricalData.from_list_df(sub_df), family_id))
        return cds

    def full_data(self):
        return [(CategoricalData.from_list_df(self.df), "full")]

    def get_languages(self):
        if self.df is None:
            return []
        return self.df['Language_ID'].unique()


class CLDFCognateData(NativeData):

    def __init__(self, source_path):
        self.df = None
        forms_path = os.path.join(source_path, "forms.csv")
        if not os.path.isfile(forms_path):
            return
        cognates_path = os.path.join(source_path, "cognates.csv")
        if not os.path.isfile(cognates_path):
            return
        forms_df = pd.read_csv(forms_path)
        if len(forms_df.index) == 0:
            return
        cognates_df = pd.read_csv(cognates_path)
        if len(cognates_df.index) == 0:
            return

        forms_df = drop_columns_except(forms_df, ["ID", "Language_ID", "Parameter_ID"])
        forms_df = forms_df.rename(columns={'Parameter_ID': 'Char_ID'})
        forms_df = forms_df.astype({'Language_ID':'string'})

        cognates_df = cognates_df[cognates_df.Cognate_Detection_Method == "expert"]
        cognates_df = drop_columns_except(cognates_df, ["Form_ID", "Cognateset_ID"])
        cognates_df = cognates_df.rename(columns={'Cognateset_ID': 'Value'})

        self.df = pd.merge(forms_df, cognates_df, how="outer", left_on=['ID'], right_on=['Form_ID'])
        self.df = self.df[self.df.Value == self.df.Value]
        self.df = self.df[self.df.Char_ID == self.df.Char_ID]
        self.df = self.df[self.df.Language_ID == self.df.Language_ID]
        self.df = self.df.drop("Form_ID", axis=1)
        self.df = self.df.drop("ID", axis=1)
        self.df = self.df.astype(str)


class CLDFStrcuturalData(NativeData):

    def __init__(self, source_path):
        self.df = None
        values_path = os.path.join(source_path, "values.csv")
        if not os.path.isfile(values_path):
            return
        self.df = pd.read_csv(values_path)
        if len(self.df.index) == 0:
            return
        self.df = drop_columns_except(self.df, ["Language_ID", "Parameter_ID", "Value"])
        self.df = self.df.rename(columns={'Parameter_ID': 'Char_ID'})
        self.df = self.df.astype({'Language_ID':'string'})

        indefinite_values = ["-", "?", "N/A", "X", "None"]
        for value in indefinite_values:
            self.df.replace(value, float("nan"))
        self.df = self.df[self.df.Value == self.df.Value]
        self.df = self.df[self.df.Char_ID == self.df.Char_ID]
        self.df = self.df[self.df.Language_ID == self.df.Language_ID]
        self.df = self.df.astype(str)


class CLDFHandler():

    def __init__(self, source_path):
        self.source_path = source_path


    def get_data(self, ling_type):
        native_data = self.get_native_data(ling_type)
        if native_data is None:
            return []
        return native_data.get_data(glottolog.split_families(self.get_glottocodes(native_data.get_languages())[0]))

    def get_sha(self):
        with open(os.path.join(self.source_path, "meta.json"), 'r') as openfile:
            json_data = json.load(openfile)
        return json_data["sha"]

    def get_dialect_ratio(self, ling_type):
        native_data = self.get_native_data(ling_type)
        glottocodes, missing_codes = self.get_glottocodes(native_data.get_languages())
        num_glottocodes = len(glottocodes)
        num_languages = native_data.num_taxa()
        num_languages_with_glottocode = num_languages - missing_codes
        if num_languages_with_glottocode == 0:
            return float('nan')
        assert(num_languages_with_glottocode >= num_glottocodes)
        return num_glottocodes / num_languages_with_glottocode


    def get_glottocodes(self, languages_with_data):
        if len(languages_with_data) == 0:
            return({}, False)
        df = self.get_languages_df()
        if df is None:
            return({}, False)
        #if languages_with_data != []:
        #    df = df[df['ID'].isin(languages_with_data)]
        num_missing_codes = 0
        glottocodes = {}
        for index, row in df.iterrows():
            lang_id = row["ID"]
            if lang_id not in languages_with_data:
                continue
            glottocode = row["Glottocode"]
            if glottocode !=  glottocode:
                num_missing_codes += 1
                continue
            if glottocode in glottocodes:
                glottocodes[glottocode].append(lang_id)
            else:
                glottocodes[glottocode] = [lang_id]
        return (glottocodes, num_missing_codes)



    def get_languages_df(self):
        languages_path = os.path.join(self.source_path, "languages.csv")
        if not os.path.isfile(languages_path):
            return None
        df = pd.read_csv(languages_path)
        if len(df.index) == 0:
            return None
        df = drop_columns_except(df, ["ID", "Glottocode"])
        df = df.astype({'ID':'string'})
        return df



    def get_native_data(self, ling_type):
        if ling_type == "cognate":
            return CLDFCognateData(self.source_path)
        elif ling_type == "structural":
            return CLDFStrcuturalData(self.source_path)



    def format_part(self, part):
        if part.startswith("["):
            part = part[1:]
            part = part.split("]")[0]
        return part

    def readme_to_dict(self):
        readme_info = {}
        path = os.path.join(self.source_path, "README.md")
        if not os.path.isfile(path):
            return readme_info
        lines = open(path, "r").readlines()
        i = 0
        while (i < len(lines)):
            line = lines[i]
            if len(line) > 0  and line[0] == "[":
                break
            i += 1
        while (i < len(lines)):
            line = lines[i]
            if len(line) == 0  or line[0] != "[":
                break
            parts = line.split(" | ")
            if not len(parts) == 2:
                print(line)
            identifier = self.format_part(parts[0])
            info = self.format_part(parts[1])
            readme_info[identifier] = info
            i += 1
        return readme_info

    def get_source(self):
        d = self.readme_to_dict()
        if not "dc:bibliographicCitation" in d:
            return ""
        return d["dc:bibliographicCitation"]
