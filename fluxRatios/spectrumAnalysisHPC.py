import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from responseCorrectedSpecrtum import * 
from datetime import datetime, timedelta
import time
import os
from sklearn.metrics import r2_score
from addClassSpectra import class_add_l1_files_time
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import pdb

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


# def generate_element_report(out_file, elements, detected_elements, count_rates):
#     """
#     Generates a space-separated file with element name, counts per second,
#     and the ratio of counts per second with respect to Silicon (Si).

#     Parameters:
#     - out_file: Path to the output text file.
#     - elements: Dictionary of elements and their respective energies (keV).
#     - detected_elements: List of booleans indicating whether an element is detected.
#     - count_rates: Dictionary with element names as keys and counts per second as values.

#     Output:
#     - Space-separagted file with element name, counts per second, and ratios.
#     """
#     # Check if Silicon (Si) was detected\
#     Si_id = list(elements.keys()).index("Si")
#     # print(detected_elements)
#     # print(Si_id)
#     # if not detected_elements[]:
#     if not detected_elements[Si_id]:
#         print("Silicon (Si) not detected. Ratios cannot be calculated.") 
    
#         for cps in count_rates:
#             out_file.write("{0:.6f}\tNA\t".format(cps))

#         out_file.write('\n')
#         return

#     # Get the counts per second for Silicon
#     silicon_cps = count_rates[Si_id]

#     # Prepare data for detected elements
#     report_data = []
#     for element, detected in zip(elements.keys(), detected_elements):
#         if detected:
#             cps = count_rates[list(elements.keys()).index(element)]
#             ratio = cps / silicon_cps  # Ratio wrt Silicon
#             report_data.append((element, cps, ratio))

    
#     # Write data to output file 
#     det_elems = [sublist[0] for sublist in report_data if len(sublist) == 3]
#     for elem in elements:
#         if not(detected_elements[list(elements.keys()).index(elem)]):
#             out_file.write("NA\tNA\t")
#             continue
#         cts, ratioo = [(row[1], row[2]) for row in report_data if row[0]==elem][0]
#         out_file.write("{0:.6f}\t{1:.6f}\t".format(cts, ratioo))
#     out_file.write('\n')




def generate_element_report(filename, bore_lat, bore_long, iteration, out_file, elements, amplitudes, detections, fitted_mus, sigmas, r2_scores):
    """
    Generates a report for each element, including its amplitudes, detections, fitted means, sigmas,
    r2 scores, and the ratio of its amplitude versus the amplitude of silicon.

    Parameters:
    - out_file (file object): The output CSV file to write the report to.
    - elements (list): List of element names.
    - amplitudes (list): List of amplitudes for each element.
    - detections (list): List of detection values for each element.
    - fitted_mus (list): List of fitted means (mus) for each element.
    - sigmas (list): List of sigmas for each element.
    - r2_scores (list): List of R^2 scores for each element.
    """
    
    # Initialize the CSV writer
    writer = csv.writer(out_file)
    
    # # Write the header for each element
    # headers = []
    # for element in elements:
    #     headers.append(f"{element}_amplitude")
    #     headers.append(f"{element}_detection")
    #     headers.append(f"{element}_fitted_mu")
    #     headers.append(f"{element}_sigma")
    #     headers.append(f"{element}_r2_score")
    #     headers.append(f"{element}_amplitude_ratio_vs_Si")  # Assuming Si is one of the elements
    
    # # Write the headers to the file
    # writer.writerow(headers)
    
    # Assuming Silicon (Si) is an element in the list, find its index to calculate amplitude ratios
    if 'Si' in elements:
        # si_index = elements.index('Si')
        si_index = list(elements.keys()).index("Si")
    else:
        raise ValueError("Silicon ('Si') must be one of the elements in the list to calculate amplitude ratio.")
    print(amplitudes)
    row = [filename]
    row.append(bore_lat)
    row.append(bore_long)
    row.append((iteration + 1) * 8)
    
    # Write the data for each element
    for i in range(len(amplitudes)):  # Assuming amplitudes, detections, fitted_mus, sigmas, r2_scores are lists of lists
        # for j, element in enumerate(elements):
        # Get the amplitude ratio with Silicon (Si)
        if amplitudes[i] != 0 and amplitudes[si_index] != 0:
            amplitude_ratio = amplitudes[i] / amplitudes[si_index]
        else:
            amplitude_ratio = 0.0  # If either is zero, set ratio to 0
        
        # Append the data for this element to the row
        
        row.append(amplitudes[i])
        row.append(detections[i])
        row.append(fitted_mus[i])
        row.append(sigmas[i])
        row.append(r2_scores[i])
        row.append(amplitude_ratio)
    
    # Write the row to the CSV file
    writer.writerow(row)






# Define combined Gaussian function: sum of 8 Gaussians plus baseline
def combined_gaussian(x, *params):
    """
    Combined Gaussian model: sum of 8 Gaussians centered at theoretical peaks with shifts,
    amplitudes, sigmas, plus a baseline.

    Parameters:
    - x: Energy array.
    - params: List of parameters [A1, shift1, sigma1, A2, shift2, sigma2, ..., A8, shift8, sigma8, baseline]

    Returns:
    - Combined Gaussian model evaluated at x.
    """
    baseline = params[-1]
    model = baseline
    num_elements = 8
    for i in range(num_elements):
        amplitude = params[i*3]
        shift = params[i*3 + 1]
        sigma = params[i*3 + 2]
        mu = theoretical_peaks[i] + shift
        model += amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return model

def model_without_element_i(x, popt, i):
    baseline = popt[-1]
    model = baseline
    num_elements = len(element_list)
    for j in range(num_elements):
        if j == i:
            continue
        amplitude = popt[j*3]
        shift = popt[j*3 +1]
        sigma = popt[j*3 +2]
        mu = theoretical_peaks[j] + shift
        model += amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return model


def analyze_spectrum_with_plot(isPlot, energy, counts, elements, detection_threshold=2, R2_threshold=0.85):
    """
    Analyze energy spectrum by fitting the entire spectrum with a sum of Gaussians,
    plot the spectrum with Gaussian fits, and output lists of count rates and detection statuses.

    Parameters:
    - isPlot: Boolean indicating whether to plot the spectrum.
    - energy: Array of energy values.
    - counts: Array of count rates corresponding to energy.
    - elements: Dictionary of elements and their respective theoretical energies (keV).
    - detection_threshold: Minimum amplitude threshold to detect an element.

    Returns:
    - Two lists:
      - amplitudes: List of fitted amplitudes for each element.
      - detections: List of booleans indicating detection status for each element.
    """
    
    # Initialize lists for results
    amplitudes = []
    detections = []
    fitted_mus = []
    sigmas = []
    r2_scores = []

    if isPlot:
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(energy, counts, s=2, label="Observed Spectrum", color="black", linewidth=1.5)
    
    # Initial guesses for parameters: [A1, shift1, sigma1, ..., A8, shift8, sigma8, baseline]
    p0 = []
    for peak in theoretical_peaks:
        p0.extend([max(counts), 0, 0.05])  # [Amplitude, Shift, Sigma]
    p0.append(1.0)  # Initial guess for baseline
    
    # Define bounds for parameters
    lower_bounds = []
    upper_bounds = []
    for peak in theoretical_peaks:
        lower_bounds.extend([0, -0.1, 1e-2])  # [A_min, shift_min, sigma_min]
        upper_bounds.extend([max(counts)*5, 0.1, 0.1])  # [A_max, shift_max, sigma_max]
    lower_bounds.append(0)  # Baseline_min
    upper_bounds.append(10)  # Baseline_max
    
    try:
        # Perform the combined Gaussian fit
        t1_=time.time()
        popt, pcov = curve_fit(combined_gaussian, energy, counts, p0=p0, bounds=(lower_bounds, upper_bounds), maxfev=400)
        t2_ = time.time()
        print("Fitting time: ", t2_-t1_)

        # Calculate fitted counts
        fitted_counts = combined_gaussian(energy, *popt)
        
        # Calculate R2 for the overall fit
        r2Value = r2_score(counts, fitted_counts)
        print("Overall R² value:", r2Value * 100, "%")
        
        # Extract parameters for each Gaussian
        for i, elem in enumerate(element_list):
            amplitude = popt[i*3]
            shift = popt[i*3 + 1]
            sigma = popt[i*3 + 2]
            fitted_mu = theoretical_peaks[i] + shift
            print(f"{elem}: Amplitude={amplitude:.2f}, Shift={shift:.4f} keV, Sigma={sigma:.4f} keV, Center={fitted_mu:.4f} keV")
            
            # Determine detection status
            detected = amplitude > detection_threshold
            # print(detected, amplitude, detection_threshold)
            amplitudes.append(amplitude)
            detections.append(detected)
            fitted_mus.append(fitted_mu)
            sigmas.append(sigma)

            # Plot individual Gaussian if plotting is enabled
            if isPlot:
                individual_gaussian = amplitude * np.exp(-0.5 * ((energy - fitted_mu) / sigma) ** 2)
                ax.plot(energy, individual_gaussian, linestyle="--", linewidth=1.2, label=f"Fitted {elem}")
        
        # Plot the overall fit
        if isPlot:
            ax.plot(energy, fitted_counts, color="red", label="Combined Fit", linewidth=2)
    
    except Exception as e:
        print("Combined fit failed:", e)
        # If fit fails, return zeros and False detections
        amplitudes = [0.0] * num_elements
        detections = [False] * num_elements

        return [None]*5
    
    if isPlot:
        # Beautify the plot
        ax.set_xlabel("Energy (keV)", fontsize=14)
        ax.set_ylabel("Counts", fontsize=14)
        ax.set_title("Energy Spectrum with Combined Gaussian Fits", fontsize=16, fontweight="bold")
        ax.legend(title="Legend", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        plt.show()
    
    # pdb.set_trace()
    
    for i, elem in enumerate(element_list):
        amplitude = popt[i*3]
        shift = popt[i*3 + 1]
        sigma = popt[i*3 + 2]
        fitted_mu = theoretical_peaks[i] + shift

        # Create a model that is the sum of Gaussians of all OTHER elements (except itself)
        model_wo_i = model_without_element_i(energy, popt, i)

        # Subtract model_wo_i from counts data to get residuals
        residuals = counts - model_wo_i

        # Create the individual Gaussian model for the current element
        individual_gaussian = amplitude * np.exp(-0.5 * ((energy - fitted_mu) / sigma) ** 2)

        # Select data within 3*sigma around the elemental peak
        mask = np.abs(energy - fitted_mu) <= 3 * sigma
        energy_window = energy[mask]
        residuals_window = residuals[mask]
        individual_gaussian_window = individual_gaussian[mask]

        # Calculate R² score for the residuals against the individual Gaussian
        r2_element = r2_score(residuals_window, individual_gaussian_window)
        print(f"R² score for {elem} within ±3σ: {r2_element:.4f}")

        r2_scores.append(r2_element)

        detections[i] *= (0.03< sigmas[i] <0.1) and (r2_element>R2_threshold)

        if isPlot:
            # Plot residuals against the individual Gaussian model within the window
            plt.figure(figsize=(12, 8))
            plt.scatter(energy_window, residuals_window, s=2, label=f"Residuals (Counts - Model w/o {elem})", color="black")
            plt.plot(energy_window, individual_gaussian_window, label=f"Gaussian Fit for {elem}", linewidth=2, color="blue")

            # Beautify the plot
            plt.xlabel("Energy (keV)", fontsize=14)
            plt.ylabel("Residual Counts", fontsize=14)
            plt.title(f"Residuals and Gaussian Fit for {elem}", fontsize=16, fontweight="bold")
            plt.legend(title=f"R² = {r2_element:.4f}", fontsize=12)
            plt.grid(True, linestyle="--", alpha=0.7)
            plt.tight_layout()
            plt.show()
    
    return amplitudes, detections, fitted_mus, sigmas, r2_scores


def process_spectra(
    dir_path,
    rmf_file,
    arf_file,
    bg_file,
    output_file,
    elements,
    detection_threshold,
    R2_threshold,
    start_counter=1,
    end_counter=None
):
    """
    Processes spectral files within a directory, analyzes them, and generates a report.

    Parameters:
    - dir_path (str): Path to the directory containing FITS files.
    - rmf_file (str): Path to the RMF file.
    - arf_file (str): Path to the ARF file.
    - bg_file (str): Path to the background FITS file.
    - output_file (str): Path to the output report text file.
    - elements (list): List of elements to analyze.
    - detection_threshold (float): Threshold for detection in analysis.
    - start_counter (int, optional): Starting file counter. Defaults to 1.
    - end_counter (int, optional): Ending file counter. If None, processes all files from start_counter onwards.
    """
    
    # Start timing
    t1 = time.time()
    
    # List all FITS files in the directory
    file_names = list_fits_files_in_directory(dir_path)
    
    # Open the output file
    with open(output_file, 'w+') as out_file:

        # Write header
        header = "filename,BORE_LAT,BORE_LON,Integration_Time,Fe_amplitude,Fe_detection,Fe_fitted_mu,Fe_sigma,Fe_r2_score,Fe_amplitude_ratio_vs_Si, \
Ti_amplitude,Ti_detection,Ti_fitted_mu,Ti_sigma,Ti_r2_score,Ti_amplitude_ratio_vs_Si, \
Ca_amplitude,Ca_detection,Ca_fitted_mu,Ca_sigma,Ca_r2_score,Ca_amplitude_ratio_vs_Si, \
Si_amplitude,Si_detection,Si_fitted_mu,Si_sigma,Si_r2_score,Si_amplitude_ratio_vs_Si, \
Al_amplitude,Al_detection,Al_fitted_mu,Al_sigma,Al_r2_score,Al_amplitude_ratio_vs_Si, \
Mg_amplitude,Mg_detection,Mg_fitted_mu,Mg_sigma,Mg_r2_score,Mg_amplitude_ratio_vs_Si, \
Na_amplitude,Na_detection,Na_fitted_mu,Na_sigma,Na_r2_score,Na_amplitude_ratio_vs_Si, \
O_amplitude,O_detection,O_fitted_mu,O_sigma,O_r2_score,O_amplitude_ratio_vs_Si" + '\n'
        out_file.write(header)
        
        file_counter = 1
        
        # Iterate over each FITS file
        for spec_file in file_names:

            # Check if the current file should be processed based on the counters
            if file_counter < start_counter or (end_counter is not None and file_counter > end_counter):
                file_counter += 1
                continue
            hdul = fits.open(os.path.join(dir_path, spec_file))
            primary_hdu = hdul[1]
            header = primary_hdu.header
            keys = set(header.keys())
            bore_lat = header.get('BORE_LAT')
            bore_long = header.get('BORE_LON')
            print(bore_long)
            # exit()
            print(f"\n\n ============ File Number {file_counter} : {spec_file} ===============")
            current_spectrum = os.path.join(dir_path, spec_file)
            op_dir_location = os.path.join(dir_path, "INTEGRATED")
            if not os.path.exists(op_dir_location):
                os.makedirs(op_dir_location)
            # Extract time windows from filename
            filename = os.path.basename(current_spectrum)
            end_time_str = filename[-23:-5]
            start_time_str = filename[-42:-24]
            time_format = "%Y%m%dT%H%M%S%f"
            start_time = datetime.strptime(start_time_str, time_format)
            end_time = datetime.strptime(end_time_str, time_format)

            adjust_end = True
            TIME_STEP = 8
            iteration = 0
            max_integration_duration = 40
            max_iteration = max_integration_duration / TIME_STEP
            print(start_time, end_time)
            while(True):
                # Process the spectrum
                energy, counts = process_spectrum_with_fits(
                    current_spectrum,
                    rmf_file,
                    arf_file,
                    bg_file,
                    ''
                )
                
                # Analyze the spectrum
                amplitudes, detections, fitted_mus, sigmas, r2_scores = analyze_spectrum_with_plot(
                    True,
                    energy,
                    counts,
                    elements,
                    detection_threshold,
                    R2_threshold
                )

                if not(amplitudes):
                    print("Useless file!")
                    break
                
                # Generate the report for the current element
                # generate_element_report(out_file, elements, detections, amplitudes)
                generate_element_report(spec_file, bore_lat, bore_long, iteration, out_file, elements, amplitudes, detections, fitted_mus, sigmas, r2_scores)
              
                # Print the analysis results
                print(amplitudes, detections, fitted_mus, sigmas, r2_scores)

                # Adjust time windows
                if adjust_end:
                    end_time += timedelta(seconds=TIME_STEP)
                    print(f"Adjusting end time: +{TIME_STEP} seconds")
                else:
                    start_time -= timedelta(seconds=TIME_STEP)
                    print(f"Adjusting start time: -{TIME_STEP} seconds")          
                adjust_end = not adjust_end  # Toggle between end and start time adjustment
                iteration += 1
                if iteration >= max_iteration:
                    break

                t1_ = time.time()
                integrated_file_name = class_add_l1_files_time(dir_path, start_time.strftime("%Y-%m-%dT%H:%M:%S.%f"), end_time.strftime("%Y-%m-%dT%H:%M:%S.%f"), op_dir_location)
                t2_ = time.time()
                print("Add spectra time: ", t2_-t1_)
                
                
                current_spectrum = integrated_file_name

            file_counter += 1
    
    # End timing
    t2 = time.time()
    
    print(f"Processing completed in {t2 - t1:.2f} seconds.")
    return



# Elements and their corresponding energies (in keV)
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

detection_threshold = 2
R2_threshold = 0.85
element_list = list(elements.keys())
theoretical_peaks = [elements[elem] for elem in element_list]
num_elements = len(element_list)


# SINGLE FOLDER PROCESSING (WORKING)

# Define necessary variables
dir_path = "/mnt/hdd/Projects/Lunar-Mapping/testingPipeline/INTEGRATED"
rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
bg_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/background_allevents.fits"
output_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/output_68.csv"


# Define the range of files to process (e.g., files 120 to 122)
start_counter = 350
end_counter = 122

# Call the function
process_spectra(
    dir_path=dir_path,
    rmf_file=rmf_file,
    arf_file=arf_file,
    bg_file=bg_file,
    output_file=output_file,
    elements=elements,
    detection_threshold=detection_threshold,
    R2_threshold=R2_threshold,
    # start_counter=start_counter,
    # end_counter=end_counter
)





# SINGLE CORE PROCESSING (UNTESTED)


# def process_all_flares(flare_category_dir, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold):
#     # Loop through each date directory inside the flare category directory
#     for date_dir in os.listdir(flare_category_dir):
#         date_path = os.path.join(flare_category_dir, date_dir)
        
#         # Only proceed if it's a directory and follows the expected format (e.g., 20210828)
#         if os.path.isdir(date_path) and date_dir.isdigit() and len(date_dir) == 8:
            
#             # Loop through each flare_XX folder inside the date directory
#             for flare_folder in os.listdir(date_path):
#                 flare_path = os.path.join(date_path, flare_folder)
                
#                 # Only process directories that match the pattern "flare_XX"
#                 if os.path.isdir(flare_path) and flare_folder.startswith('flare_'):
#                     flare_number = flare_folder.split('_')[1]  # Extract the flare number (XX)
#                     output_file = os.path.join(date_path, f"output_{flare_number}.csv")  # Define output file
                    
#                     # Call the process_spectra function for each flare_XX folder
#                     process_spectra(
#                         dir_path=flare_path,
#                         rmf_file=rmf_file,
#                         arf_file=arf_file,
#                         bg_file=bg_file,
#                         output_file=output_file,
#                         elements=elements,
#                         detection_threshold=detection_threshold,
#                         R2_threshold=R2_threshold
#                     )
#                     print(f"Processed {flare_path} -> Output saved to {output_file}")
#         else:
#             print(f"Skipping non-date directory: {date_dir}")

# # Example usage:
# flare_category_dir = "/home/sammy/Downloads/Processing/M"
# rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
# arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
# bg_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/background_allevents.fits"

# # Call the function to process all flares
# process_all_flares(
#     flare_category_dir=flare_category_dir,
#     rmf_file=rmf_file,
#     arf_file=arf_file,
#     bg_file=bg_file,
#     elements=elements,
#     detection_threshold=detection_threshold,
#     R2_threshold=R2_threshold
# )


# def process_flare(flare_path, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold, output_file):
#     """
#     Helper function to process a single flare directory.
#     """
#     process_spectra(
#         dir_path=flare_path,
#         rmf_file=rmf_file,
#         arf_file=arf_file,
#         bg_file=bg_file,
#         output_file=output_file,
#         elements=elements,
#         detection_threshold=detection_threshold,
#         R2_threshold=R2_threshold
#     )
#     return f"Processed {flare_path} -> Output saved to {output_file}"

# def process_all_flares_parallel(flare_category_dir, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold):
#     # Create a list to store the tasks
#     tasks = []

#     # Set up the ProcessPoolExecutor with up to 4 workers
#     with ProcessPoolExecutor(max_workers=4) as executor:
#         # Loop through each date directory inside the flare category directory
#         for date_dir in os.listdir(flare_category_dir):
#             date_path = os.path.join(flare_category_dir, date_dir)

#             # Only proceed if it's a directory and follows the expected format (e.g., 20210828)
#             if os.path.isdir(date_path) and date_dir.isdigit() and len(date_dir) == 8:
                
#                 # Loop through each flare_XX folder inside the date directory
#                 for flare_folder in os.listdir(date_path):
#                     flare_path = os.path.join(date_path, flare_folder)
                    
#                     # Only process directories that match the pattern "flare_XX"
#                     if os.path.isdir(flare_path) and flare_folder.startswith('flare_'):
#                         flare_number = flare_folder.split('_')[1]  # Extract the flare number (XX)
#                         output_file = os.path.join(date_path, f"output_{flare_number}.csv")  # Define output file
                        
#                         # Submit the task to the executor
#                         tasks.append(executor.submit(process_flare, flare_path, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold, output_file))

#         # Wait for all futures to complete and handle results
#         for future in as_completed(tasks):
#             print(future.result())

# Example usage:
# flare_category_dir = "/home/sammy/Downloads/Processing/M"
# rmf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_rmf_v1.rmf"
# arf_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/class_arf_v1.arf"
# bg_file = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/background_allevents.fits"
# script_path = "/mnt/hdd/Projects/Lunar-Mapping/fluxRatios/"
# elements = ["element1", "element2"]  # Define as needed
# detection_threshold = 0.5  # Adjust based on your needs
# R2_threshold = 0.95  # Adjust based on your needs

# Call the function to process all flares in parallel
# process_all_flares_parallel(
#     flare_category_dir=flare_category_dir,
#     rmf_file=rmf_file,
#     arf_file=arf_file,
#     bg_file=bg_file,
#     elements=elements,
#     detection_threshold=detection_threshold,
#     R2_threshold=R2_threshold
# )


# def process_flare_in_terminal(flare_path, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold, output_file, script_path):
#     """
#     This function launches a new terminal window to execute the processing of a flare.
#     """
#     # Command to run the process_spectra function in a new terminal
#     env_path = "/mnt/hdd/Projects/Lunar-Mapping/python_env"
#     command = [
#         'alacritty', '-e', 'bash', '-c',
#         f'echo "Processing {flare_path}"; '
#         f'echo $PYTHONPATH; '
#         f'export PYTHONPATH=$PYTHONPATH:{os.path.dirname(script_path)}; '
#         f'source {env_path}/bin/activate; '

#         f'which python; '
#         f'echo "Python envt selected!"; '
#         f'echo $PYTHONPATH; '

#         f'python3 -c "from spectrumAnalysis import process_spectra; '
#         f'process_spectra(\'{flare_path}\', \'{rmf_file}\', \'{arf_file}\', \'{bg_file}\', '
#         f'\'{output_file}\', {elements}, {detection_threshold}, {R2_threshold})"; '
#         'echo "Press any key to close the terminal..."; read -n 1'
#     ]
#         # Open a new terminal and execute the command
#     subprocess.Popen(command)

#         # f'python3 -c "from spectrumAnalysis import process_spectra; '
#         # f'python -c "from spectrumAnalysis import process_spectra; '

# def process_all_flares_parallel_with_terminal(flare_category_dir, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold, script_path):
#     # Loop through each date directory inside the flare category directory
#     for date_dir in os.listdir(flare_category_dir):
#         date_path = os.path.join(flare_category_dir, date_dir)

#         # Only proceed if it's a directory and follows the expected format (e.g., 20210828)
#         if os.path.isdir(date_path) and date_dir.isdigit() and len(date_dir) == 8:
            
#             # Loop through each flare_XX folder inside the date directory
#             for flare_folder in os.listdir(date_path):
#                 flare_path = os.path.join(date_path, flare_folder)
                



                
#                 # Only process directories that match the pattern "flare_XX"
#                 if os.path.isdir(flare_path) and flare_folder.startswith('flare_'):
#                     flare_number = flare_folder.split('_')[1]  # Extract the flare number (XX)
#                     output_file = os.path.join(date_path, f"output_{flare_number}.csv")  # Define output file
                    
#                     # Launch a new terminal for processing the flare
#                     # pdb.set_trace()
#                     # if input() in ['X', 'x']: exit()
#                     process_flare_in_terminal(flare_path, rmf_file, arf_file, bg_file, elements, detection_threshold, R2_threshold, output_file, script_path)
#                     print(f"Processing {flare_path} in a new terminal window...")


# if __name__ == "__main__":
#     process_all_flares_parallel_with_terminal(
#         flare_category_dir=flare_category_dir,
#         rmf_file=rmf_file,
#         arf_file=arf_file,
#         bg_file=bg_file,
#         elements=elements,
#         detection_threshold=detection_threshold,
#         R2_threshold=R2_threshold,
#         script_path=script_path
#     )
