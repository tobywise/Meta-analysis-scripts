__author__ = 'k1327409'
import nibabel
from scipy import ndimage
import numpy as np
import os
import re


def check_jk_niftis(mean_niftis, jk_dir):
    """
    Compares jack-knife nifti outputs with an original nifti file to check whether clusters disappear

    Arguments:
    ----------
    mean_niftis = List of thresholded positive and negative mean analysis niftis
    jk_dir = Directory containing the thresholded results of the jackknife analysis

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

    # structure for labeling
    s = [[[1,1,1],
         [1,1,1],
         [1,1,1]],
        [[1,1,1],
         [1,1,1],
         [1,1,1]],
        [[1,1,1],
         [1,1,1],
         [1,1,1]]]

    studies = []
    jk_niftis = []
    for jk in jk_files:
        if re.match(r'.+_z_p_.+(?<!p\.nii)\.gz', jk):
            #print jk
            jk_niftis.append(jk)
            study_name = re.findall(r'(?<=[a-z]JK).+_et_al_', jk)[0]
            studies.append(study_name)

    studies = set(studies)
    studies = list(studies)

    results_dict = {}

    for study in studies:
        jk_iter = []
        for jk in jk_niftis:
            if study in jk:
                jk_iter.append(jk)

        for j in jk_iter:
            if '_neg' in j:
                jk_nifti_n = nibabel.load(jk_dir + '/' + j).get_data()
            else:
                jk_nifti_p = nibabel.load(jk_dir + '/' + j).get_data()

        jk_nifti_p[np.where(jk_nifti_p == 0)] = -999  # sets anything with 0 value in JK to -999
        jk_nifti_n[np.where(jk_nifti_n == 0)] = -999  # sets anything with 0 value in JK to -999
        jk_nifti_p[np.where(original_nifti_p == 0)] = 0  # mask with original image - anything 0 in original becomes 0
        jk_nifti_n[np.where(original_nifti_n == 0)] = 0  # mask with original image - anything 0 in original becomes 0

        print "Checking " + study
        for img in ['p', 'n']:
            # Label masked JK image
            if img == 'p':
                jk_img = jk_nifti_p
                original_img = original_nifti_p
                original_aff = original_aff_p
            elif img == 'n':
                jk_img = jk_nifti_n
                original_img = original_nifti_n
                original_aff = original_aff_n

            labeled_array, num_features = ndimage.label(jk_img, structure=s)

            # NEED TO SOMEHOW LABEL CLUSTERS SO YOU KNOW WHICH ONE IS WHICH
            for i in np.unique(labeled_array):
                if i != 0:
                    print "Cluster " + str(i)
                    m = np.mean(jk_img[np.where(labeled_array == i)])
                    missing_vox = float(np.sum(jk_img[np.where(labeled_array == i)] == -999))
                    total_vox = float(labeled_array[np.where(labeled_array == i)].shape[0])

                    max = np.max(original_img[np.where(labeled_array == i)])
                    max_coords = np.where((original_img == max) & (labeled_array == i))
                    max_coords_arr = np.zeros(3)
                    for i in range(0, len(max_coords)):
                        max_coords_arr[i] = max_coords[i][0]
                    max_coords_mni = nibabel.affines.apply_affine(original_aff, max_coords_arr)
                    if not str(max_coords_mni) in results_dict:
                        results_dict[str(max_coords_mni)] = []
                    print max_coords_mni
                    print "Cluster size = " + str(int(total_vox)) + " voxels"
                    print "Percent missing voxels = " + str(round(missing_vox/total_vox*100, 2)) + "%"
                    if m == -999:
                        print "Cluster not present"
                        results_dict[str(max_coords_mni)].append(study)

    return results_dict

# Example
mean_niftis = ["C:/Users/k1327409/Documents/VBShare/0502_mean_BD_z_p_0.00500_1.000_10.nii.gz",
               "C:/Users/k1327409/Documents/VBShare/0502_mean_BD_z_p_0.00500_1.000_10_neg.nii.gz"]
jk_dir = 'C:/Users/k1327409/Documents/VBShare/BD JK'

jackknife_check_output = check_jk_niftis(mean_niftis, jk_dir)

