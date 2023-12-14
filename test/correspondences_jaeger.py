#%%
import pandas as pd
import numpy as np
import os
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

#%%

def bin_align(in_path, frequency_repeat = False):
    d = pd.read_table(in_path)
    languages = d.columns[3:-2]
    languages = languages.sort_values()
    new_cols = d.columns[:3].union(languages, sort = False)
    new_cols = new_cols.union(d.columns[-2:], sort = False)
    d = d[new_cols]
    d = d.sort_values("ID")

    matrices = []
    s = 0
    for c in d.index:
        char = d.iloc[c,:]
        #valueTokens = char.values[3:-2].astype('unicode')
        if char.FREQUENCY > 2 or frequency_repeat:
            valueTokens = char.values[3:-2]
            valueTypes = pd.unique(valueTokens)
            valueTypes = pd.Series([s for s in valueTypes if s != 'Ø'])
            valueTypes = valueTypes.sort_values()
            if frequency_repeat:
                cMtx = np.repeat(np.transpose(np.vstack([valueTokens == s for s in valueTypes], dtype=int)), int(char.FREQUENCY), axis=1)
                #assert(np.shape(cMtx)[1] == len(valueTypes) * int(char.FREQUENCY))
                #s += len(valueTypes) * int(char.FREQUENCY)
            else:
                cMtx = np.repeat(np.transpose(np.vstack([valueTokens == s for s in valueTypes], dtype=int)), 1, axis=1)
                #assert(np.shape(cMtx)[1] == len(valueTypes))
                #s += len(valueTypes)
            cMtx[valueTokens == 'Ø',:] = -1
            matrices.append(cMtx)
    charMtx = pd.DataFrame(np.hstack(matrices), index=languages)
    #assert(len(charMtx.columns) == s)
    characters = np.array(['0', '1', '-'])
    sequences = []
    for i, row in charMtx.iterrows():
        sequences.append("".join([characters[r] for r in row]))
    #for sequence in sequences:
        #assert(len(sequence) == s)
    records = [SeqRecord(sequences[i],id=language) for (i, language) in enumerate(languages)]
    msa = MultipleSeqAlignment(records, annotations={}, column_annotations={})
    return msa
# %%
