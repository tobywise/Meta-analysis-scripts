__author__ = 'k1327409'
import nibabel
from scipy import ndimage
import numpy as np

original_nifti = nibabel.load('C:/Users/k1327409/Documents/VBShare/0502_mean_BD_z_p_0.00500_1.000_10.nii.gz').get_data()

jk_nifti = nibabel.load('C:/Users/k1327409/Documents/VBShare/BD JK/1602_BD_Mean_JckJKcombined_t_Ha_et_al_2009_sMRI_z_p_0.00500_1.000_10.nii.gz').get_data()
jk_aff = nibabel.load('C:/Users/k1327409/Documents/VBShare/BD JK/1602_BD_Mean_JckJKcombined_t_Ha_et_al_2009_sMRI_z_p_0.00500_1.000_10.nii.gz').get_affine()

jk_nifti[np.where(original_nifti == 0)] = -99

s = [[[1,1,1],
     [1,1,1],
     [1,1,1]],
    [[1,1,1],
     [1,1,1],
     [1,1,1]],
    [[1,1,1],
     [1,1,1],
     [1,1,1]]]

# NEXT STEP IS TO WORK OUT WHAT I NEED TO LABEL, ORIGINAL OR JK
labeled_array, num_features = ndimage.label(jk_nifti, structure=s) # I THINK IT WORKS, GIVES 8 CLUSTERS WHICH IS SAME AS POS SECTION OF HTML!!!!!!

for i in np.unique(labeled_array):
    if i != 0:
        print i
        m = np.mean(jk_nifti[np.where(labeled_array == i)])
        print m
        if m == 0:
            print "Cluster not present"

new_image = nibabel.Nifti1Image(jk_nifti, jk_aff)
nibabel.save(new_image, "testimg.nii.gz")

masked_jk = mask_data(jk_nifti, original_nifti)




