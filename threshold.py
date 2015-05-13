import nibabel
import numpy as np


def threshold_img(original_image, threshold):
    """
    Thresholds a nifti image

    Arguments:
    original_image = image to be thresholded (in .nii.gz format)
    threshold = threshold (e.g. 0.005)

    Returns:
    Nothing, just saves a thresholded version

    E.g.
    threshold_img('path/to/file.nii.gz', 0.005)

    """
    data = nibabel.load(original_image).get_data()
    affine = nibabel.load(original_image).get_affine()

    data[np.where(data < (1-threshold))] = 0

    img = nibabel.Nifti1Image(data, affine)

    img.to_filename(original_image.replace('.nii.gz', '_' + str(threshold) + '.nii.gz'))

    print "Thresholded " + original_image + " at " + str(threshold)

