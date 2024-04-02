import os
import pandas as pd
import lingdata.pathbuilder as pb
import lingdata.params as params
import lingdata.segmetrics as segmetrics
from lingpy.align import pairwise

def drop_columns_except(df, relevant_columns):
    new_df = df
    for column in df.columns:
        if not column in relevant_columns:
            new_df = new_df.drop(column, axis=1)
    return new_df



def generate_membership_msa(ds_id, source, ling_type, family, dist_metric):
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
            all_segemnts_lists = []
            languages = list(class_df['Language_ID'].unique())
            languages.sort()
            for language in languages:
                language_df = class_df[class_df["Language_ID"] == language]
                segments_list = list(language_df['Segments'].unique())
                if len(segments_list) == 1 and segments_list[0] == "nan":
                    return
                segments_list = [segments.split(" ") for segments in segments_list]
                segments_list = [list(filter(lambda s: s != "+", segments)) for segments in segments_list]
                segments_list = [segmetrics.remove_slashes(segments) for segments in segments_list]
                all_segemnts_lists.append(segments_list)
            dm = [[0 for i in range(len(all_segemnts_lists))] for j in range(len(languages))]
            for i, segments_list1 in enumerate(all_segemnts_lists):
                for j in range(i):
                    segments_list2 = all_segemnts_lists[j]
                    distances = []
                    for s1 in segments_list1:
                        for s2 in segments_list2:
                            if dist_metric == "lev":
                                distances.append(segmetrics.levenshtein(s1, s2))
                            elif dist_metric == "jaro":
                                distances.append(segmetrics.jaro(s1, s2))
                            elif dist_metric == "mattis":
                                pair = pairwise.Pairwise(s1, s2)
                                pair.align(method='sca', distance=True)
                                distances.append(pair.alignments[0][-1])
                    d = sum(distances) / len(distances) #maybe also min or max?
                    #print(segments_list1)
                    #print(segments_list2)
                    #print(d)
                    dm[i][j] = d
                    dm[j][i] = d
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
                one_prob = 1 - zero_prob
                one_prob = round(one_prob, 3)
                zero_prob = round(zero_prob, 3)
                prob_vec.append(str(zero_prob)+","+str(one_prob))
            cols.append(col)
            probs.append(prob_vec)




    with open(pb.msa_path(ds_id, source, "cognate", "full", "membership_" + dist_metric), "w+") as outfile:
        outfile.write(str(len(all_languages)) + " ")
        outfile.write(str(num_sites) + "\n")
        outfile.write(" ".join(all_languages))
        outfile.write("\n")
        for i in range(num_sites):
            outfile.write(cols[i] + " ")
            outfile.write(" ".join(probs[i]))
            outfile.write("\n")
