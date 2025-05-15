# debugging
# print("I am inside responsecorrectedspectrum!")
# import subprocess
# command = [
#         'alacritty', '-e', 'bash', '-c',
#         f'which python; '
#         f'echo "inside responsecorrectedspectrum!"; '
#     ]
# subprocess.Popen(command)

from typing import clear_overloads
# from xspec import *
from astropy.io import fits
import numpy as np

def process_spectrum_with_fits(spec_file, rmf_file, arf_file, bg_fits_file, output_file):
    """
    Processes an energy spectrum file with RMF, ARF, and background FITS file.
    Outputs corrected energy vs. count rate values into a text file.

    Parameters:
    - spec_file: Path to the spectrum file (PHA).
    - rmf_file: Path to the response matrix file (RMF).
    - arf_file: Path to the auxiliary response file (ARF).
    - bg_fits_file: Path to the background file (FITS).
    - output_file: Path to the output text file.
    """
    try:
        # AllData.clear()
        # Xset.chatter = 0
        # Load the spectrum
        # spectrum = Spectrum(spec_file)

        # Assign response matrix and ARF files
        # spectrum.response = rmf_file
        # spectrum.response.arf = arf_file

        # Get the spectrum exposure time
        # spec_exposure = spectrum.exposure

        # Load spectrum data from FITS file
        with fits.open(spec_file) as hdul:
            spec_data = hdul[1].data  # Assuming background data is in the first extension
            spec_channels = spec_data['CHANNEL']
            spec_counts = spec_data['COUNTS']
            spec_exposure = hdul[1].header['EXPOSURE']  # spectrum exposure time
            
        # Load background data from FITS file
        with fits.open(bg_fits_file) as hdul:
            bg_data = hdul[1].data  # Assuming background data is in the first extension
            # bg_channels = bg_data['CHANNEL']
            bg_counts = bg_data['COUNTS']
            bg_exposure = hdul[1].header['EXPOSURE']  # Background exposure time

        # Ensure the response is loaded
        # if not spectrum.response:
        #     raise ValueError("Response files must be properly loaded.")

        # Plot data to initialize the energy and counts arrays
        # Plot.xAxis = "keV"  # Set x-axis to energy
        # Plot.background = True
        # Plot("data")  # Plot to initialize the data arrays

        # Extract energies and spectrum counts
        # energies = np.array(Plot.x())
        # spec_counts = np.array(Plot.y())

        # Match background counts to spectrum channels
        # if len(bg_channels) != len(spec_counts):
        #     raise ValueError("Background channels do not match spectrum channels.")

        # Convert counts to count rates (counts per second)
        spec_count_rate = spec_counts / spec_exposure
        bg_count_rate = bg_counts / bg_exposure

        # Perform background correction
        corrected_count_rate = spec_count_rate - bg_count_rate

        # Ensure non-negative corrected count rates
        corrected_count_rate = np.maximum(corrected_count_rate, 0)

        energies = np.array(spec_channels) * 13.5 / 1000


        return energies, corrected_count_rate
    
        # # Write corrected energy vs. count rate to output file
        # with open(output_file, "w") as f:
        #     f.write("# Energy (keV)\tCorrected Count Rate (counts/s)\n")
        #     for e, c in zip(energies, corrected_count_rate):
        #         f.write(f"{e:.3f}\t{c:.6f}\n")

        # print(f"Corrected energy vs. count rate written to {output_file}")

    except Exception as e:
        print(f"Error processing spectrum: {e}")
        return None, None

# # Example usage
# spec_file = "ch2_cla_l1_20210828T062321451_20210828T062329451.fits"
# rmf_file = "class_rmf_v1.rmf"
# arf_file = "class_arf_v1.arf"
# bg_file = "background_allevents.fits"
# output_file = "corrected_spectrum.txt"

# process_spectrum_with_fits(spec_file, rmf_file, arf_file, bg_file, output_file)
