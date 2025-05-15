import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from responseCorrectedSpecrtum import * 
import time
import os
from sklearn.metrics import r2_score


def list_files_in_directory(directory_path):
    try:
        # Get a list of all entries in the directory
        all_entries = os.listdir(directory_path)
        
        # Filter only files
        files = [entry for entry in all_entries if os.path.isfile(os.path.join(directory_path, entry))]
        
        return files
    except FileNotFoundError:
        print(f"Error: Directory '{directory_path}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


import os

def list_fits_files_in_directory(directory_path):
    # Define the allowed FITS file extensions
    fits_extensions = {'.fits', '.fit', '.fts'}
    
    try:
        # Get a list of all entries in the directory
        all_entries = os.listdir(directory_path)
        
        # Filter only files with FITS extensions
        fits_files = [
            entry for entry in all_entries
            # if os.path.isfile(os.path.join(directory_path, entry)) and
            #    os.path.splitext(entry)[1].lower() in fits_extensions
            if os.path.splitext(entry)[1].lower() in fits_extensions
        ]
        
        return fits_files
    except FileNotFoundError:
        print(f"Error: Directory '{directory_path}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

# Define Gaussian function
def gaussian(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

# Elements and their corresponding energies (in keV, example values)
elements = {
    "Fe": 6.40,
    "Ti": 4.51,
    "Ca": 3.69,
    "Si": 1.74,
    "Al": 1.49,
    "Mg": 1.25,
    "Na": 1.04,
    "O": 0.52
}

# Detection threshold
detection_threshold = 2

def generate_element_report(output_file, elements, detected_elements, count_rates):
    """
    Generates a space-separated file with element name, counts per second,
    and the ratio of counts per second with respect to Silicon (Si).

    Parameters:
    - output_file: Path to the output text file.
    - elements: Dictionary of elements and their respective energies (keV).
    - detected_elements: List of booleans indicating whether an element is detected.
    - count_rates: Dictionary with element names as keys and counts per second as values.

    Output:
    - Space-separated file with element name, counts per second, and ratios.
    """
    # Check if Silicon (Si) was detected\
    Si_id = list(elements.keys()).index("Si")

    # if not detected_elements[]:
    if not detected_elements[Si_id]:
        print("Silicon (Si) not detected. Ratios cannot be calculated.") 
    
        for cps in count_rates:
            out_file.write("{0:.6f}\tNA\t".format(cps))

        out_file.write('\n')
        return

    # Get the counts per second for Silicon
    silicon_cps = count_rates[Si_id]

    # Prepare data for detected elements
    report_data = []
    for element, detected in zip(elements.keys(), detected_elements):
        if detected:
            cps = count_rates[list(elements.keys()).index(element)]
            ratio = cps / silicon_cps  # Ratio wrt Silicon
            report_data.append((element, cps, ratio))

    
    # Write data to output file 
    det_elems = [sublist[0] for sublist in report_data if len(sublist) == 3]
    for elem in elements:
        if not(detected_elements[list(elements.keys()).index(elem)]):
            out_file.write("NA\tNA\t")
            continue
        cts, ratioo = [(row[1], row[2]) for row in report_data if row[0]==elem][0]
        out_file.write("{0:.6f}\t{1:.6f}\t".format(cts, ratioo))
    out_file.write('\n')




# Function to fit data and check for element detection
def analyze_spectrum_with_plot(isPlot, energy, counts, elements, detection_threshold=1):
    """
    Analyze energy spectrum from a text file with energy and counts, plot the spectrum with Gaussian fits,
    and output a dictionary with element symbol as key and a list of count rate and element detection status.

    Parameters:
    - file_path: Path to the input text file.
    - elements: Dictionary of elements and their respective energies (keV).
    - detection_threshold: Minimum amplitude threshold to detect an element.

    Returns:
    - A dictionary where keys are element symbols, and values are lists:
      [count rate, detection status (True/False)].
    """

    thresholdWidth = 0.15

    # Define Gaussian function
    def gaussian(x, a, mu, sigma, c):
        # return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
        return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + c

    # Initialize dictionary for results
    # results = {}
    results = [[], []]

    # Create the plot
    if isPlot: 
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(energy, counts, s=2, label="Observed Spectrum", color="black", linewidth=1.5)

    for element, energy_peak in elements.items():
        # Select a range around the peak energy for fitting
        x, y = 0, 0

        mask = (energy > energy_peak - thresholdWidth) & (energy < energy_peak + thresholdWidth)
        x = energy[mask]
        y = counts[mask]

        # Skip if no data points in the range
        # if len(x) < 3:
        #     results[element] = [0.0, False]
        #     continue

        p0= [] 

        # Initial guesses for Gaussian fit parameters
        # p0 = [max(y), energy_peak, 0.005]  # [amplitude, mean, sigma]
        p0 = [max(y), energy_peak, 0.05, 1]  # [amplitude, mean, sigma, baseline]

        # Amplitude bounds
        amp_min = 0
        amp_max = max(y) * 5  # Adjust multiplier

        # Mean bounds
        mean_min = energy_peak - thresholdWidth
        mean_max = energy_peak + thresholdWidth

        # Sigma bounds
        sigma_min = 1e-6  # Avoid zero
        sigma_max = 7.5e-2 # Adjust based on data

        # Baseline bounds
        baseline_min = 1e-6
        baseline_max = 1

        lower_bounds = [amp_min, mean_min, sigma_min, baseline_min]
        upper_bounds = [amp_max, mean_max, sigma_max, baseline_max]

        try:
            # Perform Gaussian fitting
            popt, _ = curve_fit(gaussian, x, y, p0=p0, bounds=(lower_bounds, upper_bounds))

            #r2 score 
            # r2Value = r2_score(y, gaussian(x, popt[0], popt[1], popt[2]))
            r2Value = r2_score(y, gaussian(x, popt[0], popt[1], popt[2], popt[3]))
            # print(element, r2Value, sep=' : ')
            print(element, "R2 value:", r2Value * 100, "%")

            # Extract fitted amplitude (counts per second)
            amplitude = popt[0]

            # Determine detection status based on threshold
            detected = amplitude > detection_threshold
            # print(element, amplitude, detected, detection_threshold)

            # Add to results
            # results[element] = [amplitude, detected]
            results[0].append(amplitude)
            results[1].append(detected)
            print(*popt)
            # Plot the Gaussian fit
            y_fit = None
            y_fit = gaussian(x, *popt)

            if isPlot: 
                ax.plot(x, y_fit, label=f"Fitted {element}", linestyle="--", linewidth=1.2)

        except Exception as e:
            # Fit failed; mark as not detected
            # results[element] = [0.0, False]
            results[0].append(0.0)
            results[1].append(False)

    # Beautify the plot
    if isPlot: 
        ax.set_xlabel("Energy (keV)", fontsize=14)
        ax.set_ylabel("Counts", fontsize=14)
        ax.set_title("Energy Spectrum with Fitted Peaks", fontsize=16, fontweight="bold")
        ax.legend(title="Legend", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.show()
        # input()

    # input()
    return list(results)
    # return dict(results)



def split_dictionary_values(input_dict):
    """
    Splits the values of a dictionary into two separate lists.

    Parameters:
    - input_dict: A dictionary with values as 2-element lists.

    Returns:
    - first_elements: A list of the first elements of each value.
    - second_elements: A list of the second elements of each value.
    """
    first_elements = []
    second_elements = []

    for value in input_dict.values():
        if True:  # Ensure the value has exactly 2 elements
            first_elements.append(value[0])
            second_elements.append(value[1])
        else:
            raise ValueError("All dictionary values must be 2-element lists.")

    return first_elements, second_elements


# Example usage
dirPath = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/flare_68"
# fileNames = list_files_in_directory(dirPath)
fileNames = list_fits_files_in_directory(dirPath)

# print(fileNames)
# rmf_file = "class_rmf_v1.rmf"
# arf_file = "class_arf_v1.arf"
# bg_file = "background_allevents.fits"

rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
bg_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/background_allevents.fits"


t1 = time.time()

out_file = open("20210828-flare_68_test2.txt", 'w')
for e in elements:
    out_file.write(f"{e}_rate\t{e}_Si\t")
out_file.write('\n')

fileCounter = 1

for spec_file in fileNames:

    # if fileCounter == 100: break
    print (f"\n\n ============ File Number {fileCounter} : {spec_file} ===============")
    energy, counts = process_spectrum_with_fits(dirPath+'/'+spec_file, rmf_file, arf_file, bg_file, '')
    # results = None
    # results = analyze_spectrum_with_plot(fileCounter>120, energy, counts, elements, detection_threshold)
    cts, detections = analyze_spectrum_with_plot(fileCounter>120, energy, counts, elements, detection_threshold)

    # cts, detections = split_dictionary_values(results)
    # print(cts, detections)
    # for c, d in zip(elements, detections):
    #     print(c, " : ", int(d))

    # print(results.values())
    generate_element_report(out_file, elements, detections, cts)

    # print(f"Element report written to {out_file} for the file {spec_file}")

    fileCounter += 1
    # if not(input()):
    #     exit()
    # break

out_file.close()
t2 = time.time()

print("Time: ",t1, t2, t2-t1)


