__author__ = 'k1327409'
import nibabel
from scipy import ndimage

original_nifti = nibabel.load('').get_data

jk_nifti = nibabel.load('C:/Users/k1327409/Documents/VBShare/BD JK/1602_BD_Mean_JckJKcombined_t_Ha_et_al_2009_sMRI_z_p_0.00500_1.000_10.nii.gz').get_data()
jk_aff = nibabel.load('C:/Users/k1327409/Documents/VBShare/BD JK/1602_BD_Mean_JckJKcombined_t_Ha_et_al_2009_sMRI_z_p_0.00500_1.000_10.nii.gz').get_affine()

s = [[[1,1,1],
     [1,1,1],
     [1,1,1]],
    [[1,1,1],
     [1,1,1],
     [1,1,1]],
    [[1,1,1],
     [1,1,1],
     [1,1,1]]]
labeled_array, num_features = ndimage.label(jk_nifti, structure=s) # I THINK IT WORKS, GIVES 8 CLUSTERS WHICH IS SAME AS POS SECTION OF HTML!!!!!!

new_image = nibabel.Nifti1Image(labeled_array, jk_aff)
nibabel.save(new_image, "testimg.nii.gz")

masked_jk = mask_data(jk_nifti, original_nifti)




