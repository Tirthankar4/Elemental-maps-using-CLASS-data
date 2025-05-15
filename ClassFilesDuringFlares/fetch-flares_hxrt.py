import os
from pathlib import Path
import pandas as pd
from datetime import datetime, timedelta
import shutil
import time
from tqdm import tqdm

# Define the expected CLASS file prefix
class_file_prefix = "ch2_cla_l1_"



def get_flare_category(flare_class):
  """
  Determine flare category from flare_class_with_bg.
  """
  if pd.isna(flare_class):
      return 'UNKNOWN'
  first_char = flare_class[0].upper()
  # print(first_char)
  return first_char if first_char in ['A', 'B', 'C', 'M', 'X'] else 'OTHER'

def count_existing_flares(category_dir):
  """
  Count existing flares in a category directory.
  """
  max_number = 0
  if category_dir.exists():
      for date_dir in category_dir.glob('*'):
          if date_dir.is_dir():
              for flare_dir in date_dir.glob('flare_*'):
                  try:
                      number = int(flare_dir.name.split('_')[1])
                      max_number = max(max_number, number)
                  except (IndexError, ValueError):
                      continue
  return max_number


def fetch_class_data():
  """
  Fetch CLASS data for all flares listed in a single catalog and copy files within
  the time window to the corresponding flare directories.
  """
  # Define base directories
  base_processing_dir = "Processing_HXRT"
  base_data_dir = "/mnt/EVENmoreSTUFF/Lunar-Mapping"
  catalog_file = base_data_dir + "/flare/flare_catalog.csv"

  # Load the flare catalog
  print("Loading flare catalog...")
  flare_catalog = pd.read_csv(catalog_file)
  total_flares = len(flare_catalog)
  print(f"Total flares to process: {total_flares}")

  # Dictionary to keep track of flare numbers per category
  category_counters = {}

  # Start timing
  start_time = time.time()

  # Create progress bar
  pbar = tqdm(total=total_flares, desc="Processing flares")

  # Loop over each flare in the catalog
  for index, flare in flare_catalog.iterrows():
      # date_str = flare['0']
      flare_category = get_flare_category(flare['class'])

      # Initialize or get counter for this category
      if flare_category not in category_counters:
          category_dir = Path(base_processing_dir) / flare_category
          category_counters[flare_category] = count_existing_flares(category_dir)

      # Increment flare number for this category
      category_counters[flare_category] += 1
      current_flare_number = category_counters[flare_category]

      # Calculate ETA
      elapsed_time = time.time() - start_time
      flares_processed = index + 1
      flares_remaining = total_flares - flares_processed
      eta = (elapsed_time / flares_processed) * flares_remaining if flares_processed > 0 else 0

      # Update progress bar
      pbar.set_postfix({
          #'Category': flare_category,
          'flare#': current_flare_number,
          'ETA': f'{eta:.1f}s'
      })
      pbar.update(1)
      peak_str = flare['peak']

      # Construct the path to the CLASS data directory for the date
      year = str(peak_str)[:4]
      month = str(peak_str)[5:7]
      day = str(peak_str)[8:10]
      date_str = year + month + day
      # print(date_str)
      class_data_path = Path(base_data_dir) / "CLASS Data" / "cla" / "data" / "calibrated" / year / month / day 
      # print(class_data_path)

      # Check if CLASS data directory exists for this date
      if not class_data_path.exists():
          continue

      # Create a directory for this flare
      flare_dir = Path(base_processing_dir) / flare_category / str(date_str) / f"flare_{current_flare_number}"
      flare_dir.mkdir(parents=True, exist_ok=True)


      # Convert to datetime with assumed reference date
      flare_start = datetime.strptime(flare['start'], "%Y/%m/%d %H:%M")
      flare_end = datetime.strptime(flare['end'], "%Y/%m/%d %H:%M")
      # print(flare_start, flare_end)
      # Define time window for searching CLASS data
      start_window = flare_start - timedelta(minutes=15)
      end_window = flare_end + timedelta(minutes=15)

      # Search for CLASS data files within this time window
      for class_file in class_data_path.iterdir():
          if class_file.is_file() and class_file.name.startswith(class_file_prefix) and class_file.name.endswith(".fits"):
              try:
                  parts = class_file.stem.split('_')
                  file_start_str = parts[3]
                  file_end_str = parts[4]

                  # Convttttttttert file times to datetime
                  file_start = datetime.strptime(file_start_str, "%Y%m%dT%H%M%S%f")
                  file_end = datetime.strptime(file_end_str, "%Y%m%dT%H%M%S%f")

                  # Check if the file's time range overlaps with the flare time window
                  if (file_start <= end_window) and (file_end >= start_window):
                      # Copy the file to the flare-specific directory
                      # print("Found", class_file)
                      destination = flare_dir / class_file.name
                      # print(class_file, destination)
                      shutil.copy2(class_file, destination)
                  # else:
                  #   print("U RETARD")

              except (IndexError, ValueError) as e:
                  continue
          # else:
          #       print("U RETARD")

  pbar.close()
  print("\nProcessing Summary:")
  for category, count in category_counters.items():
      print(f"Category {category}: {count} flares processed")

def main():
  print("Starting CLASS data fetching process...")
  fetch_class_data()
  print("CLASS data fetching process completed.")

if __name__ == "__main__":
  main()