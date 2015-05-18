
# script to run an entire meta-analysis

import os
import subprocess
import shlex
import pandas as pd
import numpy as np
from sdm_functions import threshold_jackknife, get_coords, extract_coordinate_values
from check_jk_niftis import check_jk_niftis

ma_dir = 'C:/MA/thing/'
sdm_path = 'C:/sdm'
analysis_name = '18_05_MDD'
filter_var = 'CombinedGroups'
sdm_table = 'C:\Users\k1327409\Documents\VBShare/2602_analysis\Extract\sdm_table.txt'


os.chdir(ma_dir)

#  preprocess
print "Preprocessing"
subprocess.call(shlex.split(sdm_path + ' pp gray_matter, 1.0, 20, gray_matter, 20'))
#  mean
print "Running mean analysis"
subprocess.call(shlex.split(sdm_path + ' ' + analysis_name + '_mean' + ' = mean ' + filter_var))

#  threshold mean & heterogeneity
print "Thresholding mean and heterogeneity"
subprocess.call(shlex.split(sdm_path + ' threshold ' + analysis_name + '_mean_z' + ', p, 0.005, 1, 10'))
subprocess.call(shlex.split(sdm_path + ' threshold ' + analysis_name + '_mean_QH_z' + ', p, 0.005, 1, 10'))

#  jack-knife - no way to call from command line so have to do manually :(
print "Running jack-knife"
table = pd.read_csv(sdm_table, delimiter='\t')

num_studies = len(table)

selection_var = table['CombinedGroups']

for i in range(1,num_studies):
    if selection_var[i] == 1:  # only go through this if study is selected
        thing = np.ones(num_studies)  # generate an array of ones the length of the table
        thing[i] = 0  # set the current study to 0 (we want to leave this study out of this iteration)
        name = analysis_name + '_JACKKNIFE_' + table['study'][thing == 0].values[0]  # name of the analysis
        thing[np.where(selection_var == 0)] = 0  # set anything where selection variable is zero to zero in thing
        table['JK_column'] = thing  # add this array to the table in a new column
        print name  # print stuff because I like printing stuff
        print thing
        table.to_csv('sdm_table.txt', sep='\t', index=False)  # save the table with this new column
        subprocess.call([sdm_path, name, ' = ', 'mean', 'JK_column'])  # run a mean analysis using this JK column

#  run and threshold meta-regressions
regression_vars = ['Antidepressants', 'Duration', 'HAMD17', 'Patient_age', 'PercentMale']

for j in regression_vars:
    metareg_name = analysis_name + '_' + j
    subprocess.call(shlex.split(sdm_path + ' ' + metareg_name + ' = ' + j + ', ' + filter_var))
    print shlex.split(sdm_path + ' ' + metareg_name + ' = ' + j + ', ' + filter_var)
    subprocess.call(shlex.split(sdm_path + ' threshold ' + metareg_name + '_1m0_z' + ', p, 0.0005, 1, 10'))

#  threshold JKs
threshold_jackknife(ma_dir, sdm_path)

#  extract peaks from mean and meta-regressions
mean_coords = get_coords(ma_dir + analysis_name + '_mean_z_p_0.005_1.000_10.htm')
extract_coordinate_values(mean_coords, analysis_name + ' _mean', sdm_path)

for j in regression_vars:
    metareg_name = analysis_name + '_' + j
    coords = get_coords(ma_dir + metareg_name + '_1m0_z_p_0.0005_1.000_10.htm')
    extract_coordinate_values(coords, metareg_name, sdm_path)

#  check jack-knife and meta-regressions
mean_niftis = [ma_dir + analysis_name + '_mean_z_p_0.00500_1.000_10.nii.gz',
               ma_dir + analysis_name + '_mean_z_p_0.00500_1.000_10_neg.nii.gz']
qh_niftis = [ma_dir + analysis_name + '_mean_QH_z_p_0.00500_1.000_10.nii.gz',
             ma_dir + analysis_name + '_mean_QH_z_p_0.00500_1.000_10.nii.gz']

jackknife_check_output = check_jk_niftis(mean_niftis, ma_dir, csv_name=analysis_name + '_JK_check.csv')

metareg_mean_check_output = check_jk_niftis(mean_niftis, ma_dir, regex=r'^.+(?=_1m0)',
                                            csv_name=analysis_name + '_mean_metareg_check.csv', metareg=True)
metareg_qh_check_output = check_jk_niftis(qh_niftis, ma_dir, regex=r'^.+(?=_1m0)',
                                          csv_name=analysis_name + '_QH_metareg_check.csv', metareg=True)
