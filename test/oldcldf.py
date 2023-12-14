import os
import pandas as pd
import math
import numpy as np
from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import json

states = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+",
",", "/", ":", ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]

column_threshold = 0
drop_threshold = 0.1

class OldCLDFHandler:
    source_path = ""

    def __init__(self, source_path):
        self.source_path =  source_path

    def drop_unnecessary_columns(self, df):
        relevant_columns = ["ID", "Form_ID", "Language_ID", "Cognateset_ID", "Name", "Parameter_ID"]
        new_df = df
        for column in df.columns:
            if not column in relevant_columns:
                new_df = new_df.drop(column, axis=1)
        return new_df

    def get_df(self, type):
        if type == "cognate":
            df = self.df_for_cognate()
        elif type == "structural":
            df = self.df_for_structural()
        else:
            raise ValueError("Type " + type + " not defined")
        return df


    def get_bin_align(self, ling_type):
        df = self.get_df(ling_type)
        languages = list(set(df['Language_ID']))
        matrix = {}
        for language in languages:
            matrix[language] = ""
        chars = list(set(df['Char_ID']))
        chars.sort()
        num_chars = len(chars)
        for (c, char) in enumerate(chars):
            sub_df = df[df["Char_ID"] == char]
            char_languages = list(set(sub_df['Language_ID']))
            values = list(set(sub_df['Value']))
            values.sort()
            for language in languages:
                if language not in char_languages:
                    bs = ['-' for _ in range(len(values))]
                    matrix[language] += ''.join(bs)
                    continue
                bs = ['0' for _ in range(len(values))]
                for(i, value) in enumerate(values):
                    sub_sub_df = sub_df[(sub_df["Language_ID"] == language) & (sub_df["Value"] == value)]
                    if(len(sub_sub_df.index) > 0):
                        bs[i] = '1'
                matrix[language] += ''.join(bs)
        records = [SeqRecord(matrix[language_id], id=str(language_id)) for language_id in matrix.keys()]
        align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
        return align


    def get_multi_align(self, ling_type):
        df = self.get_df(ling_type)
        languages = list(set(df['Language_ID']))
        languages.sort()
        matrix = {}
        for language in languages:
            matrix[language] = ""
        chars = list(set(df['Char_ID']))
        chars.sort()
        r = {
            "num_languages" : len(languages),
            "num_chars" : len(chars),
            "multi_cells" : 0,
            "multi_chars" : 0,
            "high_multi_chars" : 0,
            "high_state_chars" : 0,
            "dropped_chars" : 0,
            "max_states" : 0,
            "converted" : 0}
        for (c, char) in enumerate(chars):
            char_values = {}
            possible = True
            sub_df = df[df["Char_ID"] == char]
            char_languages = list(set(sub_df['Language_ID']))
            classes = list(set(sub_df['Value']))
            if (len(classes) > len(states)):
                r["high_state_chars"] += 1
                r["dropped_chars"] += 1
                continue
            counts = sub_df['Value'].value_counts()
            classes.sort()
            column_multi_cells = 0
            for language in languages:
                sub_sub_df = sub_df[(sub_df["Language_ID"] == language)]
                sub_sub_df = sub_sub_df.drop_duplicates()
                if language not in char_languages or len(sub_sub_df) == 0:
                    char_values[language] = '-'
                    continue
                max_count = 0
                if len(sub_sub_df) > 1:
                    column_multi_cells += 1
                for i, row in sub_sub_df.iterrows():
                    if(counts[row['Value']] > max_count):
                        char_values[language] =  states[classes.index(row["Value"])]
                assert(language in char_values)
            if column_multi_cells / len(languages) > column_threshold:
                r["high_multi_chars"] += 1
                r["dropped_chars"] += 1
            else:
                for language in languages:
                    matrix[language] += char_values[language]
            r["max_states"] = max(r["max_states"], len(classes))
            r["multi_cells"] += column_multi_cells
            if column_multi_cells > 0:
                r["multi_chars"] += 1
        if (r["dropped_chars"]/r["num_chars"] > drop_threshold):
            return (r, None)
        r["converted"] = 1
        records = [SeqRecord(matrix[language_id], id=str(language_id)) for language_id in matrix.keys()]
        align = MultipleSeqAlignment(records, annotations={}, column_annotations={})
        return (r, align)



    def write_catg(self, ling_type, path):
        df = self.get_df(ling_type)
        languages = list(set(df['Language_ID']))
        languages.sort()
        matrix = {}
        chars = list(set(df['Char_ID']))
        chars.sort()
        tempfile = open("temp.txt", "w+")
        tempfile.write(" ".join(languages))
        tempfile.write("\n")
        num_sites = 0
        for (c, char) in enumerate(chars):
            value_counts = {}
            for language in languages:
                value_counts[language] = 0
            char_df = df[df["Char_ID"] == char]
            values = list(set(char_df['Value']))
            values.sort()
            #determine number of values of this characteristic per language
            for(i, value) in enumerate(values):
                value_df = char_df[char_df["Value"] == value]
                value_languages = list(set(value_df['Language_ID']))
                for language in value_languages:
                    value_counts[language] += 1
            #determine most probable cols and probablities and write to file
            num_sites += len(values)
            for(i, value) in enumerate(values):
                value_df = char_df[char_df["Value"] == value]
                value_languages = list(set(value_df['Language_ID']))
                col = ""
                probs = []
                for language in languages:
                    if language in value_languages:
                        # only in this case 1 is more probable than 0
                        if value_counts[language] == 1:
                            col += "1"
                        else:
                            col += "0"
                        one_prob = 1 / value_counts[language]
                    else:
                        col += "0"
                        # missing data
                        if value_counts[language] == 0:
                            probs.append("1.0,1.0")
                            continue
                            #alternative:
                            #one_prob = 1 / len(values)
                        else:
                            one_prob = 0.0
                    zero_prob = 1.0 - one_prob
                    one_prob = round(one_prob, 3)
                    zero_prob = round(zero_prob, 3)
                    probs.append(str(zero_prob)+","+str(one_prob))
                tempfile.write(col + " ")
                tempfile.write(" ".join(probs))
                tempfile.write("\n")
        outfile = open(path, "w+")
        outfile.write(str(len(languages)) + " ")
        outfile.write(str(num_sites) + "\n")
        tempfile.close()
        outfile.write(open("temp.txt", "r").read())
        os.remove("temp.txt")



    def df_for_cognate(self):
        forms_path = os.path.join(self.source_path, "forms.csv")
        cognates_path = os.path.join(self.source_path, "cognates.csv")
        if not os.path.isfile(forms_path) or not os.path.isfile(cognates_path):
            return None
        forms_df = pd.read_csv(forms_path, low_memory=False)
        cognates_df = pd.read_csv(cognates_path, low_memory=False)
        cognates_df = cognates_df[cognates_df.Cognate_Detection_Method == "expert"]
        forms_df = self.drop_unnecessary_columns(forms_df)
        cognates_df = self.drop_unnecessary_columns(cognates_df)

        full_df = pd.merge(forms_df, cognates_df, how="outer", left_on=['ID'], right_on=['Form_ID'])
        full_df = full_df[full_df.Cognateset_ID == full_df.Cognateset_ID]
        full_df = full_df.rename(columns={'Cognateset_ID': 'Value', 'Parameter_ID': 'Char_ID'})
        full_df = full_df.astype(str)
        full_df = full_df.drop("ID_x", axis=1)
        full_df = full_df.drop("ID_y", axis=1)
        full_df = full_df.drop("Form_ID", axis=1)
        return full_df


    def df_for_structural(self):
        values_path =  os.path.join(self.source_path, "values.csv")
        if not os.path.isfile(values_path):
            return None
        df = pd.read_csv(values_path, low_memory=False)
        columns = ["Language_ID", "Parameter_ID", "Value"]
        new_df = df
        for column in df.columns:
            if not column in columns:
                new_df = new_df.drop(column, axis=1)
        df = new_df

        indefinite_values = ["-", "?", "N/A", "X", "None"]
        for value in indefinite_values:
            df.replace(value, np.NaN)
        df = df[df.Value == df.Value]
        df = df.rename(columns={'Parameter_ID': 'Char_ID'})
        df = df.astype(str)
        return df
