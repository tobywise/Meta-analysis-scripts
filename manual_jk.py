__author__ = 'k1327409'

import numpy as np
import pandas as pd
import subprocess
import os

######################################################################
# Very poorly written script to manually run SDM jack-knife analyses #
######################################################################

sdm_table = 'C:\Users\k1327409\Documents\VBShare/2602_analysis\Extract\sdm_table.txt'

table = pd.read_csv(sdm_table, delimiter='\t')

num_studies = len(table)

selection_var = table['CombinedGroups']

selection_var2 = np.ones(len(selection_var))
selection_var2[10:] = 0
selection_var2[0:3] = 0
os.chdir('C:/Users/k1327409/Documents/VBShare/2602_analysis/Extract')

sdm_path = 'C:/Users/k1327409/Dropbox/PhD/Things/sdm_v4.12/sdm_v4.12/sdm.bat'

for i in range(1,num_studies):
    if selection_var[i] == 1 and selection_var2[i] == 1:
        thing = np.ones(num_studies)
        thing[i] = 0
        name = 'MDD_JackJK' + table['study'][thing == 0].values[0]
        thing[np.where(selection_var == 0)] = 0  # set anything where selection variable is zero to zero in thing
        table['JK_column'] = thing
        print name
        print thing
        table.to_csv('sdm_table.txt', sep='\t', index=False)
        subprocess.call([sdm_path, name, ' = ', 'mean', 'JK_column'])
