import os
import pandas as pd


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

def load_and_process_data(element):
    # Import the CSV file
    file_path = f"Subsets/{element}_subset.csv"
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} was not found.")

    data = pd.read_csv(file_path)

    # Ensure the relevant headers exist
    required_columns = ['BORE_LAT', 'BORE_LONG', f'{element}_amplitude_ratio_vs_Si']
    if not all(col in data.columns for col in required_columns):
        raise ValueError(f"Missing one or more required columns: {required_columns}")

    return data

# define a grid of 18000x36000 pixels
subpixel_array = np.zeros((18000, 36000))
