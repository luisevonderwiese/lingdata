import os
import pandas as pd
from Levenshtein import distance
import lingdata.pathbuilder as pb
import lingdata.params as params

def drop_columns_except(df, relevant_columns):
    new_df = df
    for column in df.columns:
        if not column in relevant_columns:
            new_df = new_df.drop(column, axis=1)
    return new_df



def generate_membership_msa(ds_id, source, ling_type, family):
    if ling_type != "cognate":
        return  #makes only sence for cognate data
    if family != "full":
        return #family splitting not supported
    if source not in params.source_types["cldf"]:
        return #only supported fopr cldf
    source_path = pb.source_path("native", ds_id, source)
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

    forms_df = drop_columns_except(forms_df, ["ID", "Language_ID", "Parameter_ID", "Segments"])
    forms_df = forms_df.rename(columns={'Parameter_ID': 'Char_ID'})
    forms_df = forms_df.astype({'Language_ID':'string'})

    cognates_df = cognates_df[cognates_df.Cognate_Detection_Method == "expert"]
    cognates_df = drop_columns_except(cognates_df, ["Form_ID", "Cognateset_ID"])
    cognates_df = cognates_df.rename(columns={'Cognateset_ID': 'Value'})

    forms_df = pd.merge(forms_df, cognates_df, how="outer", left_on=['ID'], right_on=['Form_ID'])
    forms_df = forms_df[forms_df.Value == forms_df.Value]
    forms_df = forms_df[forms_df.Char_ID == forms_df.Char_ID]
    forms_df = forms_df[forms_df.Language_ID == forms_df.Language_ID]
    forms_df = forms_df.drop("Form_ID", axis=1)
    forms_df = forms_df.drop("ID", axis=1)
    forms_df = forms_df.astype(str)
    forms_df.drop_duplicates(keep=False, inplace=True)
    concepts = list(forms_df['Char_ID'].unique())
    all_languages = list(forms_df['Language_ID'].unique())
    all_languages.sort()
    num_sites = 0
    cols = []
    probs = []
    for concept in concepts:
        concept_df = forms_df[forms_df["Char_ID"] == concept]
        cognate_classes = list(concept_df['Value'].unique())
        num_sites += len(cognate_classes)
        for cognate_class in cognate_classes:
            class_df = concept_df[concept_df['Value'] == cognate_class]
            forms = []
            languages = list(class_df['Language_ID'].unique())
            languages.sort()
            for language in languages:
                language_df = class_df[class_df["Language_ID"] == language]
                language_forms = [segments.replace(" ", "") for segments in list(language_df['Segments'].unique())]
                forms.append(language_forms)
            dm = [[0 for i in range(len(languages))] for j in range(len(languages))]
            for i, forms_lang1 in enumerate(forms):
                for j in range(i):
                    forms_language2 = forms[j]
                    distances = []
                    for form1 in forms_lang1:
                        for form2 in forms_language2:
                            l = max(len(form1), len(form2))
                            distances.append(distance(form1, form2) / l)
                    d = sum(distances) / len(distances) #maybe also min or max?
                    dm[i][j] = d
                    dm[j][i] = d
            i = 0
            col = ""
            prob_vec = []
            for language in all_languages:
                if language not in languages:
                    col += "0"
                    zero_prob = 1.0
                else:
                    zero_prob = sum(dm[i]) / len(dm[i])
                    if zero_prob > 0.5:
                        col+= "0"
                    else:
                        col+= "1"
                    i += 1
                one_prob = 1 - zero_prob
                one_prob = round(one_prob, 3)
                zero_prob = round(zero_prob, 3)
                prob_vec.append(str(zero_prob)+","+str(one_prob))
            cols.append(col)
            probs.append(prob_vec)




    with open(pb.msa_path(ds_id, source, "cognate", "full", "membership"), "w+") as outfile:
        outfile.write(str(len(all_languages)) + " ")
        outfile.write(str(num_sites) + "\n")
        outfile.write(" ".join(all_languages))
        outfile.write("\n")
        for i in range(num_sites):
            outfile.write(cols[i] + " ")
            outfile.write(" ".join(probs[i]))
            outfile.write("\n")
