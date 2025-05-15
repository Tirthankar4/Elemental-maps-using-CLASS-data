import os
import re
import csv
from datetime import datetime, timedelta

def parse_timestamp(timestamp_str):
    """
    Parses the timestamp from the filename and returns a datetime object.
    The timestamp format is assumed to be YYYYMMDDThhmmssSSS, where SSS are milliseconds.
    """
    # Regular expression to extract date and time components
    match = re.match(r"(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})(\d{3})", timestamp_str)
    if not match:
        raise ValueError(f"Invalid timestamp format: {timestamp_str}")
    year, month, day, hour, minute, second, millisecond = match.groups()
    dt = datetime(
        int(year), int(month), int(day),
        int(hour), int(minute), int(second),
        int(millisecond) * 1000  # Convert milliseconds to microseconds
    )
    return dt

def read_night_times(csv_file):
    """
    Reads the night_time.csv file and returns a list of night intervals.
    Each night interval is a tuple of (start_time, end_time).
    """
    night_intervals = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            start_time = datetime.strptime(row['20_min_before'], '%Y-%m-%d %H:%M:%S')
            end_time = datetime.strptime(row['20_min_after'], '%Y-%m-%d %H:%M:%S')
            night_intervals.append((start_time, end_time))
        
        #write night_intervals for debugging purposes:
        # with open("night_interval.csv", 'w', newline='') as csvfile:
        #     writer = csv.writer(csvfile)
        #     writer.writerow(['night intervals'])
        #     for file in night_intervals:
        #         writer.writerow([file])
        # print(f"Night intervals have been written to night_interval.csv")
        #print(night_intervals)
    return night_intervals

def fetch_night_data_files(class_file_path, bg_int_time_hours, night_times_csv, data_root):
    """
    Fetches CLASS data files from the night intervals within a time window around the given CLASS file.
    """

    # Extract CLASS filename
    class_filename = os.path.basename(class_file_path)

    # Extract the first timestamp from the CLASS filename
    match = re.search(r"_(\d{8}T\d{9})_", class_filename)
    if not match:
        raise ValueError(f"Invalid CLASS filename format: {class_filename}")
    class_timestamp_str = match.group(1)

    # Parse the CLASS file timestamp
    class_time = parse_timestamp(class_timestamp_str)

    # Calculate the time window
    half_window = timedelta(hours=bg_int_time_hours) / 2
    print(half_window)
    start_window = class_time - half_window
    print(start_window)
    end_window = class_time + half_window
    print(end_window)

    # Read night times from CSV
    night_intervals = read_night_times(night_times_csv)

    # Filter night intervals within the time window
    selected_night_intervals = []
    for start_time, end_time in night_intervals:
        # Check if the night interval falls within the time window
        if (start_time >= start_window and start_time <= end_window) or \
           (end_time >= start_window and end_time <= end_window) or \
           (start_time <= start_window and end_time >= end_window):
            selected_night_intervals.append((start_time, end_time))
    print(f"selected night intervals: {selected_night_intervals}")
    no_intervals=len(selected_night_intervals)
    #print(f"no of intervals: {no_intervals}")
    tit=no_intervals*40
    print(f"total integration time: {tit} minutes") 

    # List to store found data files
    data_files = []

    # Loop over selected night intervals
    for interval_start, interval_end in selected_night_intervals:

        # Generate list of dates within the interval
        dates_in_interval = get_dates_in_window(interval_start.date(), interval_end.date())
        for date in dates_in_interval:
            date_path = os.path.join(data_root, date.strftime('%Y'), date.strftime('%m'), date.strftime('%d'))
            if not os.path.isdir(date_path):
                continue  # Skip if directory does not exist

            # List all .fits files in the directory
            for filename in os.listdir(date_path):
                if not filename.endswith('.fits'):
                    continue

                # Extract the timestamp from the filename
                file_match = re.search(r"_(\d{8}T\d{9})_", filename)
                if not file_match:
                    continue

                file_timestamp_str = file_match.group(1)
                file_start_time = parse_timestamp(file_timestamp_str)

                # Check if the file's start time is within the night interval
                if interval_start <= file_start_time <= interval_end:
                    data_files.append(os.path.join(date_path, filename))
    return data_files, class_time

def get_dates_in_window(start_date, end_date):
    """
    Returns a list of dates between start_date and end_date inclusive.
    """
    delta = end_date - start_date
    return [start_date + timedelta(days=i) for i in range(delta.days + 1)]

# Example usage
if __name__ == "__main__":
    classDir = 'CLASS Data/cla/data/calibrated'

    # has to be done for all the class files corresponding to all the flares
    class_file = 'ch2_cla_l1_20191009T114350598_20191009T114358598.fits' 
    bg_int_time = 4  # in hours
    night_times_csv = 'night_times_cleaned.csv'  # Path to the CSV file containing night times

    #night_data_files contains the path of all CLASS files that are in the night interval
    night_data_files, class_time = fetch_night_data_files(class_file, bg_int_time, night_times_csv, classDir)
    
    # Generate the CSV filename
    interval_str = f"{bg_int_time}h"
    class_time_str = class_time.strftime('%Y%m%dT%H%M%S')
    csv_filename = f"night_{interval_str}_{class_time_str}.csv"

    # Write the list of files to the CSV
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['File_Path'])
        for file in night_data_files:
            writer.writerow([file])
    print(f"Night data files have been written to {csv_filename}")
