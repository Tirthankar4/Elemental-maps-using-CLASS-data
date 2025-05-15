# Import necessary modules
import pyspextools.io
import matplotlib.pyplot as plt
import sys
import pdb

def correct_spectrum(source_file, background_file, arf_file, rmf_file, output_file=None):
    """
    Corrects an X-ray spectrum using ARF, RMF, and background files.

    Parameters:
    - source_file (str): Path to the source spectrum file (FITS format).
    - background_file (str): Path to the background spectrum file (FITS format).
    - arf_file (str): Path to the ARF file.
    - rmf_file (str): Path to the RMF file.
    - output_file (str, optional): Path to save the corrected spectrum (FITS format).

    Returns:
    - corrected_spectrum (pyspextools.Spectrum): The corrected Spectrum object.
    """

    pdb.set_trace()

    # Load the source spectrum
    print("Loading source spectrum...")
    sourceSpectrum = pyspextools.io.Pha()
    sourceSpectrum.read(source_file)

    # Load the background spectrum
    print("Loading background spectrum...")
    background = pyspextools.io.Pha()
    background.read(background_file)

    # Subtract the background from the source
    print("Subtracting background...")
    sourceSpectrum.data -= background.data
    print("Background subtracted.")

    # Load the ARF and RMF files
    print("Loading ARF file...")
    arf = pyspextools.io.Arf().read(arf_file)
    print("Loading RMF file...")
    rmf = pyspextools.io.Rmf().read(rmf_file)

    # Apply the response files to correct the spectrum
    print("Applying ARF and RMF corrections...")
    corrected_spectrum = pyspextools.spectra.apply_response(sourceSpectrum, rmf, arf)
    print("Corrections applied.")

    # Save the corrected spectrum if output_file is specified
    if output_file:
        print(f"Saving corrected spectrum to {output_file}...")
        pyspextools.io.write_pha(corrected_spectrum, output_file)
        print("Corrected spectrum saved.")

    # Plot the corrected spectrum
    print("Plotting the corrected spectrum...")
    plt.figure(figsize=(10, 6))
    plt.step(corrected_spectrum.energy, corrected_spectrum.flux, where='mid', label='Corrected Spectrum')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Flux (erg/cm²/s/keV)')
    plt.title('Corrected X-ray Spectrum')
    plt.legend()
    plt.grid(True)
    plt.show()

    return corrected_spectrum

if __name__ == "__main__":
    # Define file paths (replace with your actual file paths)
    source_spectrum = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/ch2_cla_l1_20210828T062257451_20210828T063257451.fits"
    background_spectrum = "/mnt/hdd/Projects/Lunar-Mapping/Archive/background_allevents.fits"
    arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
    rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
    output_corrected = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/corrected_spectrum.fits"

    # # Check if all necessary files are provided
    # print(sys.argv)
    # if len(sys.argv) < 5:
    #     print("Usage: python correct_spectrum.py <source_spectrum.fits> <background_spectrum.fits> <response.arf> <response.rmf> [output_corrected.fits]")
    #     sys.exit(1)

    # source_spectrum = sys.argv[1]
    # background_spectrum = sys.argv[2]
    # arf_file = sys.argv[3]
    # rmf_file = sys.argv[4]
    # output_corrected = sys.argv[5] if len(sys.argv) > 5 else None

    # Correct the spectrum
    corrected = correct_spectrum(source_spectrum, background_spectrum, arf_file, rmf_file, output_corrected)
