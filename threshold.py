import nibabel
import numpy as np
from scipy import ndimage, stats


def threshold_img(p_image, z_image, p_threshold, z_threshold, extent):
    """
    Thresholds a nifti image

    Arguments:
    p_image = raw p value image to be thresholded (in .nii.gz format)
    z_image = z value image to be thresholded (in .nii.gz format)
    p_threshold = p value threshold (e.g. 0.005)
    z_threshold = z value threshold (e.g. 0.638)
    extent = extent threshold in voxels (e.g. 10)

    Returns:
    Nothing, just saves a thresholded version

    E.g.
    p = 'path\to\file\MDD_pos_BD_pos_praw.nii.gz'
    z = 'path\to\file\MDD_pos_BD_pos_z.nii.gz'

    threshold_img(p, z, 0.005, 0.0638, 10)

    """
    p_data = nibabel.load(p_image).get_data()
    z_data = nibabel.load(z_image).get_data()
    affine = nibabel.load(p_image).get_affine()

    p_data[np.where(p_data < (1-p_threshold))] = 0

    z_data[np.where(p_data == 0)] = 0
    z_data[np.where(z_data < z_threshold)] = 0

    # structure for labeling - there's probably a better way to define this
    s = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]

    labeled_array, num_features = ndimage.label(z_data, structure=s)  # Label clusters in z image

    sizes = stats.itemfreq(labeled_array)

    for i in range(1, np.max(labeled_array)):
        if sizes[i][1] < extent:
            z_data[np.where(labeled_array == i)] = 0

    z_img = nibabel.Nifti1Image(z_data, affine)
    z_filename = z_image.replace('.nii.gz', '_thresholded_%s_%s_%s.nii.gz' % (p_threshold, z_threshold, extent))
    z_img.to_filename(z_filename)

    print "Thresholded %s at %s, %s, %s" % (z_image, p_threshold, z_threshold, extent)

