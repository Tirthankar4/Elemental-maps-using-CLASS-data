import numpy as np
from astropy.io import fits
from astropy.convolution import convolve
from scipy.interpolate import interp1d
import pdb 

# Load the FITS file (spectrum)
def load_fits_data(fits_file):
    with fits.open(fits_file) as hdul:
        data = hdul[1].data
        energy = data['CHANNEL']  # Energy values (in keV, for example)
        counts = data['COUNTS']  # Counts per bin
    return energy, counts

# Load ARF file (Auxiliary Response Function)
def load_arf_file(arf_file):
    with fits.open(arf_file) as hdul:
        arf_data = hdul[1].data
        arf_energy = arf_data['ENERG_LO']  # ARF energy bins
        arf_values = arf_data['SPECRESP']  # ARF values at each energy bin
    return arf_energy, arf_values

# Load RMF file (Redistribution Matrix File)
def load_rmf_file(rmf_file):
    with fits.open(rmf_file) as hdul:
        rmf_data = hdul[2].data
        pdb.set_trace()
        rmf_matrix = rmf_data['MATRIX']  # Redistribution matrix
        print(rmf_matrix.shape)
        rmf_energy = rmf_data['ENERG_LO']  # RMF energy bins
    return rmf_energy, rmf_matrix

# Load Background file (if applicable)
def load_background_file(bkg_file):
    with fits.open(bkg_file) as hdul:
        bkg_data = hdul[1].data
        bkg_counts = bkg_data['COUNTS']  # Background counts
    return bkg_counts

# Correct the spectrum using ARF, RMF, and Background
def correct_spectrum(energy, counts, arf_energy, arf_values, rmf_energy, rmf_matrix, bkg_counts):
    # Interpolate ARF values for the spectrum energy bins
    arf_interp = interp1d(arf_energy, arf_values, kind='linear', bounds_error=False, fill_value=0.0)
    arf_interp_values = arf_interp(energy)
    
    # Apply RMF correction (this can be more complex depending on the structure of RMF data)
    # Assuming RMF matrix is a 2D array and we need to apply it to the spectrum
    corrected_counts = np.dot(rmf_matrix, counts)
    
    # Subtract background counts from the corrected counts
    corrected_counts -= bkg_counts
    
    # Normalize by ARF
    corrected_counts /= arf_interp_values

    return corrected_counts, energy

# Main process
def main(fits_file, arf_file, rmf_file, bkg_file):
    # Load data
    energy, counts = load_fits_data(fits_file)
    arf_energy, arf_values = load_arf_file(arf_file)
    rmf_energy, rmf_matrix = load_rmf_file(rmf_file)
    bkg_counts = load_background_file(bkg_file)

    # Correct the spectrum
    corrected_counts, corrected_energy = correct_spectrum(energy, counts, arf_energy, arf_values, rmf_energy, rmf_matrix, bkg_counts)

    # Print the corrected counts and energies
    print("Corrected Energy:", corrected_energy)
    print("Corrected Counts per Second:", corrected_counts)

# Example file paths (adjust to your specific files)
fits_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/ch2_cla_l1_20210828T062257451_20210828T063257451.fits"
arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
bkg_file = "/mnt/hdd/Projects/Lunar-Mapping/Archive/background_allevents.fits"

# Run the main function
main(fits_file, arf_file, rmf_file, bkg_file)

