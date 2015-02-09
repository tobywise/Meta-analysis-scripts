# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
import re
import subprocess
import os

# ############################################################################################################
# ## This script takes the SDM cluster results in HTML format and extracts values for the peak coordinates  ##
# ## Function 1 (get_coords) produces a list of coordinates from the HTML file                              ##
# ## Function 2 (extract_coordinate_values) takes a list of coordinates and extracts values at these points ##
# ## SDM will output the results in the current working directory - make sure to set this before running    ##
# ## Extracted files are numbered - need to make a note of which number corresponds to each cluster         ##
# ############################################################################################################


def get_coords(html_file):
    """
    Takes SDM's html results file and gets the cluster peak coordinates

    Arguments
    ---------
    html_file: The html file to read

    Returns
    -------
    coords_list: List of extracted coordinates
    """

    output = BeautifulSoup(open(html_file, 'r'))  # opens the html and makes it look nice
    # use regular expressions to find the main cluster coordinates
    coords = re.findall(r'/table>.*?-?\d+,-?\d+,-?\d+', str(output))
    # all main coordinates seem to come after a </table> bit

    coords_list = []
    for i in range(len(coords)):
        match = re.search(r'-?\d+,-?\d+,-?\d+', coords[i])  # get the actual coordinates out of the text
        coords_list.append(match.group())  # add these coordinates to a new list
    return coords_list


def extract_coordinate_values(coordinates, prefix, sdm_path):
    """
    Takes coordinates and extracts values from them using SDM

    Arguments
    ---------
    coordinates: A list of coordinates to be extracted
    prefix: A prefix for the output files (e.g. Bipolar_mean_analysis)
    sdm_path: Path to SDM (must be to the SDM.bat file)

    Returns
    -------
    Creates SDM masks and extracts coordinates for these masks
    SDM's outputs from these functions are saved
    Prints each cluster name and coordinates when run
    """

    for i in range(len(coordinates)):
        name = prefix + "_coords_" + str(i+1)  # format = prefix_coords_#
        mask_arg = name + " = mask coordinate " + coordinates[i]  # arguments for masking function
        extract_arg = "extract " + name  # arguments for extraction function
        print mask_arg  # prints mask arguments to make sure they're correct
        subprocess.call([sdm_path, mask_arg])
        subprocess.call([sdm_path, extract_arg])


# Example usage
"""
result_coords = get_coords('C:/Users/Tobyyy/Dropbox/PhD/DTI/dti_ma/MDD Only/14_03_14_MDD_Meds_1m0_z_p0.00050_1.000_10.htm')
extract_coordinate_values(result_coords, "BD_mean", "C:/Users/k1327409/Dropbox/PhD/Things/sdm_v4.12/sdm_v4.12/sdm.bat")

"""

# ############################################################################################################
# ## This function checks jackknife output files for any results that differ from the original analysis     ##
# ## ...It's a bit of a mess                                                                                ##
# ############################################################################################################


def check_jackknife(mean_results, jk_directory, save_results=True):
    """
    Checks jackknife results against original results and reports any changes

    Arguments
    ---------
    mean_results: The original HTML results file
    jk_directory: Directory where the jackknife results HTML files are located
    save_results: Save outputs to a text file (optional)

    Returns
    -------
    Prints/writes to text original coordinates, each jackknife iteration and
    coordinates that differ (new additions, original ones missing)
    """

    files = os.listdir(jk_directory)
    real_coords = get_coords(mean_results)
    results = []
    for file in files:
        if re.match(r'.+JK.+_z_.+.htm', file):
            results.append(jk_directory + '/' + file)
    print real_coords
    print results

    if save_results:
        out = open("Jackknife_check.txt", 'w')

    for result in results:
        print result
        jk_coords = get_coords(result)
        print jk_coords
        if save_results:
            out.write("Original coordinates = " + str(real_coords) + "\n")
            out.write(result + " - coordinates = " + str(jk_coords) + "\n")
        for coord in jk_coords:
            if coord not in real_coords:
                print coord + " = new cluster"
                if save_results:
                    out.write(coord + " = new cluster" + "\n")
        for real_coord in real_coords:
            if real_coord not in jk_coords:
                print real_coord + " missing from this iteration"
                if save_results:
                    out.write(real_coord + " missing from this iteration" + "\n")

# Example usage
"""
check_jackknife("J:/MA test stuff/MDD Only/14_03_14_MDD_Mean_z_p0.00500_1.000_10.htm", 'J:/MA test stuff/MDD Only')

"""
