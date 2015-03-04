import nibabel
from scipy import ndimage
import numpy as np
import os
import re


def check_jk_niftis(mean_niftis, jk_dir, regex=r'(?<=[a-z,_]JK).+(?=_z_p)'):
    """
    Compares jack-knife nifti outputs with an original nifti file to check whether clusters disappear

    Arguments:
    ----------
    mean_niftis = List of thresholded positive and negative mean analysis niftis - i.e. ['pos.nii.gz', 'neg.nii.gz']
    jk_dir = Directory containing the thresholded results of the jackknife analysis
    regex (Optional) = regex expression for identifying jack-knife outputs, can be changed to identify different outputs

    Returns:
    --------
    results_dict = A dictionary containing the peak coordinate of each cluster and studies where it isn't present

    Also prints out the results for each study including the percentage of voxels that are missing

    """

    jk_files = os.listdir(jk_dir)

    original_nifti_p = nibabel.load(mean_niftis[0]).get_data()
    original_nifti_n = nibabel.load(mean_niftis[1]).get_data()
    original_aff_p = nibabel.load(mean_niftis[0]).get_affine()
    original_aff_n = nibabel.load(mean_niftis[1]).get_affine()

    # structure for labeling - there's probably a better way to define this
    s = [[[1, 1, 1],
          [1, 1, 1],
          [1, 1, 1]],
         [[1, 1, 1],
          [1, 1, 1],
          [1, 1, 1]],
         [[1, 1, 1],
          [1, 1, 1],
          [1, 1, 1]]]

    studies = []
    jk_niftis = []
    for jk in jk_files:
        if re.match(r'.+_z_p_.+(?<!p\.nii)\.gz', jk):  # only thresholded niftis
            jk_niftis.append(jk)
            print jk
            study_name = re.search(regex, jk)  # get study name out of file name
            if study_name:
                studies.append(study_name.group())

    studies = set(studies)  # make list values unique
    studies = list(studies)

    results_dict = {}

    for study in studies:
        jk_iter = []
        for jk in jk_niftis:  # look for jk files that match the study name
            if study in jk:
                jk_iter.append(jk)

        for j in jk_iter:  # assign files to variables
            if '_neg' in j:
                jk_nifti_n = nibabel.load(jk_dir + '/' + j).get_data()
            else:
                jk_nifti_p = nibabel.load(jk_dir + '/' + j).get_data()

        jk_nifti_p[np.where(jk_nifti_p == 0)] = -999  # sets anything with 0 value in JK to -999
        jk_nifti_n[np.where(jk_nifti_n == 0)] = -999  # sets anything with 0 value in JK to -999
        jk_nifti_p[np.where(original_nifti_p == 0)] = 0  # mask with original image - anything 0 in original becomes 0
        jk_nifti_n[np.where(original_nifti_n == 0)] = 0  # mask with original image - anything 0 in original becomes 0

        print "*****************\nChecking " + study + '\n*****************'
        for img in ['p', 'n']:  # do this for both positive and negative results

            if img == 'p':
                jk_img = jk_nifti_p
                original_img = original_nifti_p
                original_aff = original_aff_p
                print "Positive Clusters\n*****************"
            else:
                jk_img = jk_nifti_n
                original_img = original_nifti_n
                original_aff = original_aff_n
                print "*****************\nNegative Clusters\n*****************"

            labeled_array, num_features = ndimage.label(jk_img, structure=s)  # Label clusters in masked JK image


            for i in np.unique(labeled_array):  # iterate over clusters
                if i != 0:
                    print "Cluster " + str(i)
                    m = np.mean(jk_img[np.where(labeled_array == i)])  # calculate mean of the cluster in the jk img
                    missing_vox = float(np.sum(jk_img[np.where(labeled_array == i)] == -999))  # number of missing vox
                    total_vox = float(labeled_array[np.where(labeled_array == i)].shape[0])

                    max = np.max(original_img[np.where(labeled_array == i)])  # get maximum value in cluster
                    max_coords = np.where((original_img == max) & (labeled_array == i))  # find where this max is
                    max_coords_arr = np.zeros(3)  # np.where outputs as tuple of arrays- need to convert to array
                    for i in range(0, len(max_coords)):
                        max_coords_arr[i] = max_coords[i][0]
                    max_coords_mni = nibabel.affines.apply_affine(original_aff, max_coords_arr)  # convert to MNI coords
                    if not str(max_coords_mni) in results_dict:
                        results_dict[str(max_coords_mni)] = []  # create empty dictionary entry for coordinates
                    print max_coords_mni
                    print "Cluster size = " + str(int(total_vox)) + " voxels"
                    print "Percent missing voxels = " + str(round(missing_vox/total_vox*100, 2)) + "%"
                    if m == -999:  # if mean = -999, indicates that the cluster is all missing
                        print "Cluster not present"
                        results_dict[str(max_coords_mni)].append(study)  # add study name to coordinate entry

    return results_dict

# Example
"""
mean_niftis = ["C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_z_p_0.00500_1.000_10.nii.gz",
               "C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_z_p_0.00500_1.000_10_neg.nii.gz"]
jk_dir = 'C:/Users/k1327409/Documents/VBShare/2602_analysis/JK/Thresholded'

jackknife_check_output = check_jk_niftis(mean_niftis, jk_dir)
"""

# Or for checking overlap of meta-regressions

mean_niftis = ["C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_z_p_0.00500_1.000_10.nii.gz",
               "C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_z_p_0.00500_1.000_10_neg.nii.gz"]
qh_niftis = ["C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_QH_z_p_0.00500_1.000_10.nii.gz",
               "C:/Users/k1327409/Documents/VBShare/2602_analysis/2602_BD_mean_QH_z_p_0.00500_1.000_10.nii.gz"]
metareg_dir = 'C:/Users/k1327409/Documents/VBShare/2602_analysis'

metareg_mean_check_output = check_jk_niftis(mean_niftis, metareg_dir, regex=r'^.+(?=_1m0)')
metareg_qh_check_output = check_jk_niftis(qh_niftis, metareg_dir, regex=r'^.+(?=_1m0)')