from __future__ import division
import nibabel
import numpy as np
from numpy import genfromtxt
import pandas as pd
from math import gamma, sqrt
import re

# CHANGE THIS TO YOUR INPUT CSV VILE
studies = pd.read_csv('C:\Users\k1327409\Google Drive\PhD\MDD BD Meta-analysis\SDM MA Files\Combined groups\combined_studies_mdd.csv')


def to_t(d, n1, n2):
    n = n1 + n2
    k = gamma((n-3) / 2) / gamma((n-2) / 2) * sqrt((n-2) * n1 * n2 / (2*n))
    t = d * k
    return t


for i in studies['study']:  # this is so messy
    print 'Study ' + str(i+1) + ' of ' + str(len(studies['study']))
    print 'Map 1 = ' + studies['group_a'][i]
    print 'Map 2 = ' + studies['group_b'][i]
    if studies['nc'][i] > 0:
        print 'Map 3 = ' + studies['group_c'][i]
    print 'N1 = ' + str(studies['na'][i]) + ', N2 = ' + str(studies['nb'][i]) + ', N3 = ' + str(studies['nc'][i])

    # CHANGE THIS PATH TO WHEREVER YOUR MAPS ARE
    maps_path = 'C:\Users\k1327409\Google Drive\PhD\MDD BD Meta-analysis\SDM MA Files\MDD sMRI/test3/'
    group1_map = maps_path + 'pp_' + studies['group_a'][i] + '.nii.gz'
    group2_map = maps_path + 'pp_' + studies['group_b'][i] + '.nii.gz'
    if studies['nc'][i] > 0:
        group3_map = maps_path + 'pp_' + studies['group_c'][i] + '.nii.gz'
    n1 = float(studies['na'][i])
    n2 = float(studies['nb'][i])
    n3 = 0
    control_n = float(studies['c_n'][i])
    if studies['nc'][i] > 0:
        n3 = float(studies['nc'][i])
    n = n1 + n2 + n3
    print n

    # load group1 nifti/affine/header
    group1 = nibabel.load(group1_map)
    group1_data = group1.get_data()
    group1_aff = group1.get_affine()
    group1_hdr = group1.get_header()

    # load group2 nifti/affine/header

    group2 = nibabel.load(group2_map)
    group2_data = group2.get_data()
    group2_aff = group2.get_affine()

    # load group 3 if it exists
    if studies['nc'][i] > 0:
        group3 = nibabel.load(group3_map)
        group3_data = group3.get_data()
        group3_aff = group3.get_affine()

    # do magic stats stuff

    # create combined effect size map
    combined_data = (n1/n) * group1_data + (n2/n) * group2_data

    if studies['nc'][i] > 0:
        combined_data = (n1/n) * group1_data + (n2/n) * group2_data + (n3/n) * group3_data

    print combined_data

    # create t-map
    combined_data = to_t(combined_data, n, control_n)

    # save combined t-map

    new_hdr = group1_hdr
    new_hdr.set_intent('t test', (n-2,))

    img = nibabel.Nifti1Image(combined_data, group1_aff, header=new_hdr)
    # CHANGE THE BIT BELOW TO WHEREVER YOU WANT THE OUTPUT TO GO
    out_path = 'C:\Users\k1327409\Google Drive\PhD\MDD BD Meta-analysis\SDM MA Files\Combined groups/'
    filename = studies['group_a'][i]
    filename = 'combined_t_' + filename.replace(re.findall('\d+.', filename)[0], re.findall('\d+', filename)[0]) + '_test.nii.gz'
    img.to_filename(out_path + filename)
