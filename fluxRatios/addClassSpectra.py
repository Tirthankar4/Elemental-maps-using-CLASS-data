import os
import glob
import numpy as np
from astropy.io import fits
from astropy.time import Time
from datetime import datetime

def convert_time_str_format(time_str):
    """Convert the file time string to the expected format."""
    yr1 = time_str[0:4]
    mon1 = time_str[4:6]
    date1 = time_str[6:8]
    hr1 = time_str[9:11]
    min1 = time_str[11:13]
    sec1 = time_str[13:15]
    mse1 = time_str[15:18]
    UTC_all = f"{yr1}-{mon1}-{date1}T{hr1}:{min1}:{sec1}.{mse1}"
    return UTC_all

def utc_to_seconds(utc):
    """Convert UTC time string to seconds since Jan 01, 1970."""
    if '.' in utc:
        time_part, ms_part = utc.split('.')
        ms_part = (ms_part + '000000')[:6]  # Pad milliseconds to microseconds
        utc_full = time_part + '.' + ms_part
        dt = datetime.strptime(utc_full, "%Y-%m-%dT%H:%M:%S.%f")
    else:
        dt = datetime.strptime(utc, "%Y-%m-%dT%H:%M:%S")
    timestamp = dt.timestamp()
    return timestamp

def utc_to_timestr(utc):
    """Convert UTC string to filename time string format."""
    timestr = utc.replace('-', '').replace(':', '').replace('.', '')
    return timestr








def write_summed_fits(summed_SCD_opfname, data_SCD_summed, ip_file, expotime, starttime, endtime, SPICE_val):
    """Write FITS file for the specified summed data."""
    # Create primary HDU
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['DATE'] = Time.now().isot  # File creation date

    # Create columns for binary table
    channels = np.arange(len(data_SCD_summed))  # Should be 2048 elements
    col1 = fits.Column(name='CHANNEL', format='J', array=channels)
    col2 = fits.Column(name='COUNTS', format='E', array=data_SCD_summed)
    cols = fits.ColDefs([col1, col2])

    # Create binary table HDU
    bintable_hdu = fits.BinTableHDU.from_columns(cols)
    bthdr = bintable_hdu.header

    # Compute MID_UTC as the average of STARTIME and ENDTIME
    start_dt = datetime.strptime(starttime[0], "%Y-%m-%dT%H:%M:%S.%f")
    end_dt = datetime.strptime(endtime[0], "%Y-%m-%dT%H:%M:%S.%f")
    mid_dt = start_dt + (end_dt - start_dt) / 2
    mid_utc = mid_dt.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]  # Format to match original

    # Add header parameters
    bthdr.update({
        'TTYPE1':    ('CHANNEL', 'PHA channel'),
        'TTYPE2':    ('COUNTS', 'Counts per channel'),
        'TUNIT2':    ('count', 'Physical unit of field'),
        'EXTNAME':   ('SPECTRUM', 'Name of binary table extension'),
        'HDUCLASS':  ('OGIP', 'format conforms to OGIP standard'),
        'HDUCLAS1':  ('SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
        'HDUVERS1':  ('1.1.0', 'Obsolete - included for backwards compatibility'),
        'HDUVERS':   ('1.1.0', 'Version of format (OGIP memo OGIP-92-007a)'),
        'HDUCLAS2':  ('UNKNOWN', 'Maybe TOTAL, NET or BKG Spectrum'),
        'HDUCLAS3':  ('COUNT', 'PHA data stored as Counts (not count/s)'),
        'TLMIN1':    (0, 'Lowest legal channel number'),
        'TLMAX1':    (2047, 'Highest legal channel number'),
        'TELESCOP':  ('CHANDRAYAAN-2', 'mission/satellite name'),
        'INSTRUME':  ('CLASS', 'instrument/detector name'),
        'FILTER':    ('none', 'filter in use'),
        'EXPOSURE':  (expotime, 'exposure (in seconds)'),
        'AREASCAL':  (1.0, 'area scaling factor'),
        'BACKFILE':  ('NONE', 'associated background filename'),
        'BACKSCAL':  (1.0, 'background file scaling factor'),
        'CORRFILE':  ('NONE', 'associated correction filename'),
        'CORRSCAL':  (1.0, 'correction file scaling factor'),
        'RESPFILE':  ('', 'associated redistrib matrix filename'),
        'ANCRFILE':  ('', 'associated ancillary response filename'),
        'PHAVERSN':  ('1992a', 'obsolete'),
        'DETCHANS':  (2048, 'total number possible channels'),
        'CHANTYPE':  ('PHA', 'channel type (PHA, PI etc)'),
        'POISSERR':  (True, 'Poissonian errors to be assumed'),
        'STAT_ERR':  (0, 'no statistical error specified'),
        'SYS_ERR':   (0, 'no systematic error specified'),
        'GROUPING':  (0, 'no grouping of the data has been defined'),
        'QUALITY':   (0, 'no data quality information specified'),
        'EQUINOX':   (2000.0, 'Equinox of Celestial coord system'),
        'DATE':      (Time.now().isot, 'file creation date'),
        'PROGRAM':   ('addClassSpectra.py', 'Program that created the file'),
        'SW_VERSN':  (2.1, 'L0-L1 processing software version'),
        'IPFILE':    (ip_file, 'Input file name'),
        'STARTIME':  (starttime[0], 'Start time in UTC'),
        'ENDTIME':   (endtime[0], 'End time in UTC'),
        'TEMP':      (round(SPICE_val['temperature_mean']*10.0) / 10.0, 'Temperature mean (degrees) rounded off to 1 dec - Mean value for the added files'),
        'GAIN':      (round(SPICE_val['gain_mean']*10000.0) / 10000.0, 'eV/channel - Mean value for the added files'),
        'SCD_FLTR':  (round(SPICE_val['scd_filter_mean']), '1 if device filtering done; 0 if not done - Mean value for the added files'),
        'SCD_USED':  ('0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,', 'SCDs used to add'),
        'MID_UTC':   (mid_utc, 'Mid window UTC'),
        'LST_HR':    (SPICE_val['LST_HR'], 'Average Local Solar Time Hour'),
        'LST_MIN':   (SPICE_val['LST_MIN'], 'Average Local Solar Time Minute'),
        'LST_SEC':   (SPICE_val['LST_SEC'], 'Average Local Solar Time Second'),
        'SAT_ALT':   (round(SPICE_val['sat_alt_mean'] * 10000.0) / 10000.0, 'Sub-satellite point altitude (km) - Mean value for the added files'),
        'SAT_LAT':   (round(SPICE_val['sat_lat_mean'] * 10000.0) / 10000.0, 'Sub-satellite point latitude (deg) - Mean value for the added files'),
        'SAT_LON':   (round(SPICE_val['sat_lon_mean'] * 10000.0) / 10000.0, 'Sub-satellite point longitude (deg) - Mean value for the added files'),
        'BORE_LAT':  (round(SPICE_val['bore_lat_mean'] * 10000.0) / 10000.0, 'Boresight latitude (deg) - Mean value for the added files'),
        'BORE_LON':  (round(SPICE_val['bore_lon_mean'] * 10000.0) / 10000.0, 'Boresight longitude (deg) - Mean value for the added files'),
        'V0_LAT':    (round(SPICE_val['v0_lat'] * 10000.0) / 10000.0, 'Pixel corner 0 latitude (deg)'),
        'V1_LAT':    (round(SPICE_val['v1_lat'] * 10000.0) / 10000.0, 'Pixel corner 1 latitude (deg)'),
        'V2_LAT':    (round(SPICE_val['v2_lat'] * 10000.0) / 10000.0, 'Pixel corner 2 latitude (deg)'),
        'V3_LAT':    (round(SPICE_val['v3_lat'] * 10000.0) / 10000.0, 'Pixel corner 3 latitude (deg)'),
        'V0_LON':    (round(SPICE_val['v0_lon'] * 10000.0) / 10000.0, 'Pixel corner 0 longitude (deg)'),
        'V1_LON':    (round(SPICE_val['v1_lon'] * 10000.0) / 10000.0, 'Pixel corner 1 longitude (deg)'),
        'V2_LON':    (round(SPICE_val['v2_lon'] * 10000.0) / 10000.0, 'Pixel corner 2 longitude (deg)'),
        'V3_LON':    (round(SPICE_val['v3_lon'] * 10000.0) / 10000.0, 'Pixel corner 3 longitude (deg)'),
        'SOLARANG':  (round(SPICE_val['sol_ang_mean'] * 10000.0) / 10000.0, 'Angle between surface normal,sun vector (deg) - Mean value for the added files'),
        'PHASEANG':  (round(SPICE_val['phs_ang_mean'] * 10000.0) / 10000.0, 'Angle between boresight vector,sun vector (deg) - Mean value for the added files'),
        'EMISNANG':  (SPICE_val['emi_ang_mean'], 'Angle btwn boresight vector,emitted X-rays(deg) - Mean value for the added files'),
    })

    # Write to FITS file
    hdulist = fits.HDUList([primary_hdu, bintable_hdu])
    hdulist.writeto(summed_SCD_opfname, overwrite=True)







def class_add_l1_files_time(ipfile_dir, UTC_start_user, UTC_end_user, op_dir_location):
    """Main program to add specified L1 files based on user-provided start and end times."""
    if ipfile_dir == '' or UTC_start_user == '' or UTC_end_user == '' or op_dir_location == '':
        print("Please enter all the input values.")
        print("Usage:")
        print("CLASS_ADD_L1_FILES_TIME('/path/to/directory/having/input/L1_FITS_files/', 'yyyy-mm-ddTHH:MM:SS.MSE', 'yyyy-mm-ddTHH:MM:SS.MSE', '/path/to/directory/to/store/output/added_L1_FITS_files/')")
        return

    print('')
    print('USER INPUTS')
    print('-----------')
    print('Input files directory:', ipfile_dir)
    print('Start UTC:', UTC_start_user)
    print('End UTC:', UTC_end_user)
    print('Location to store output files:', op_dir_location)
    print('')

    # Ensure directories end with path separator
    ipfile_dir = os.path.join(ipfile_dir, '')
    op_dir_location = os.path.join(op_dir_location, '')

    fname_ext = '.fits'
    data_SCD_summed = np.zeros(2048, dtype=float)
    ip_file = ''

    # Check if output subdirectory exists, create if not
    op_subdir_name = 'L1_ADDED_FILES_TIME'
    op_subdir_path = os.path.join(op_dir_location, op_subdir_name)

    if not os.path.exists(op_subdir_path):
        print('Subdirectory (L1_ADDED_FILES_TIME) to store output files does not exist. Will create now.')
        os.makedirs(op_subdir_path)
    else:
        print('Subdirectory (L1_ADDED_FILES_TIME) to store output files exists.')

    print('')

    # Convert the input start and end UTC into seconds
    UTC_start_user_secs = utc_to_seconds(UTC_start_user)
    UTC_end_user_secs = utc_to_seconds(UTC_end_user)
    print('')

    # Read all the L1 files in the specified input directory and get the file times
    file_list = sorted(glob.glob(os.path.join(ipfile_dir, "*.fits")))
    n_files_in_ipdir = len(file_list)
    file_list_basename = [os.path.basename(f) for f in file_list]
    file_time_only = [fname[11:48] for fname in file_list_basename]
    file_start_times = [s[0:18] for s in file_time_only]
    file_end_times = [s[19:37] for s in file_time_only]

    # Convert into format needed by utc_to_seconds
    startUTC_arr = [convert_time_str_format(s) for s in file_start_times]
    endUTC_arr = [convert_time_str_format(s) for s in file_end_times]

    # Now get the start and end UTC from file in seconds
    start_UTC_arr_secs = np.array([utc_to_seconds(s) for s in startUTC_arr])
    end_UTC_arr_secs = np.array([utc_to_seconds(s) for s in endUTC_arr])

    # Find the files matching the user-provided time interval
    ind_UTC_match = np.where((start_UTC_arr_secs >= UTC_start_user_secs) & (end_UTC_arr_secs <= UTC_end_user_secs))[0]
    count_ind_UTC_match = len(ind_UTC_match)

    if count_ind_UTC_match > 0:
        # Get the first index from start and last index from end
        start_index = ind_UTC_match[0]
        end_index = ind_UTC_match[-1]

        # Get the corresponding matching UTC values
        startUTC_match_arr = [startUTC_arr[i] for i in ind_UTC_match]
        endUTC_match_arr = [endUTC_arr[i] for i in ind_UTC_match]

        UTC_start_actual = startUTC_match_arr[0]
        UTC_end_actual = endUTC_match_arr[-1]

        # Prepare to read and sum data
        files_to_add_arr = []
        sat_alt_all = []
        sol_ang_all = []
        phs_ang_all = []
        emi_ang_all = []


        temperature_all = []
        gain_all = []
        scd_filter_all = []
        sat_lat_all = []
        sat_lon_all = []
        bore_lat_all = []
        bore_lon_all = []
        lst_hr_all = []
        lst_min_all = []
        lst_sec_all = []



        # Corner coordinates from first and last files
        v0_lat1 = v1_lat1 = v2_lat1 = v3_lat1 = 0.0
        v0_lon1 = v1_lon1 = v2_lon1 = v3_lon1 = 0.0
        v0_lat2 = v1_lat2 = v2_lat2 = v3_lat2 = 0.0
        v0_lon2 = v1_lon2 = v2_lon2 = v3_lon2 = 0.0

        for idx, j in enumerate(ind_UTC_match):
            start_UTC_str = utc_to_timestr(startUTC_arr[j])
            end_UTC_str = utc_to_timestr(endUTC_arr[j])
            filename = 'ch2_cla_l1_' + start_UTC_str + '_' + end_UTC_str + fname_ext
            files_to_add_arr.append(filename)
            
            filepath = os.path.join(ipfile_dir, filename)
            if not os.path.exists(filepath):
                print('File not found:', filepath)
                continue
            with fits.open(filepath) as hdulist:
                hdr = hdulist[1].header  # Assuming data is in extension 1
                data = hdulist[1].data
                counts = data['COUNTS']
                data_SCD_summed += counts
                ip_file = hdr.get('IPFILE', '')
                sat_alt_all.append(hdr.get('SAT_ALT', 0.0))
                sol_ang_all.append(hdr.get('SOLARANG', 0.0))
                phs_ang_all.append(hdr.get('PHASEANG', 0.0))
                emi_ang_all.append(hdr.get('EMISNANG', 0.0))

                temperature_all.append(hdr.get('TEMP', 0.0))
                gain_all.append(hdr.get('GAIN', 0.0))
                scd_filter_all.append(hdr.get('SCD_FLTR', 0.0))

                sat_lat_all.append(hdr.get('SAT_LAT', 0.0))
                sat_lon_all.append(hdr.get('SAT_LON', 0.0))
                bore_lat_all.append(hdr.get('BORE_LAT', 0.0))
                bore_lon_all.append(hdr.get('BORE_LON', 0.0))

                lst_hr_all.append(hdr.get('LST_HR', 0))
                lst_min_all.append(hdr.get('LST_MIN', 0))
                lst_sec_all.append(hdr.get('LST_SEC', 0))




                if idx == 0:
                    # First file
                    v0_lat1 = hdr.get('V0_LAT', 0.0)
                    v1_lat1 = hdr.get('V1_LAT', 0.0)
                    v2_lat1 = hdr.get('V2_LAT', 0.0)
                    v3_lat1 = hdr.get('V3_LAT', 0.0)
                    v0_lon1 = hdr.get('V0_LON', 0.0)
                    v1_lon1 = hdr.get('V1_LON', 0.0)
                    v2_lon1 = hdr.get('V2_LON', 0.0)
                    v3_lon1 = hdr.get('V3_LON', 0.0)
                if idx == count_ind_UTC_match - 1:
                    # Last file
                    v0_lat2 = hdr.get('V0_LAT', 0.0)
                    v1_lat2 = hdr.get('V1_LAT', 0.0)
                    v2_lat2 = hdr.get('V2_LAT', 0.0)
                    v3_lat2 = hdr.get('V3_LAT', 0.0)
                    v0_lon2 = hdr.get('V0_LON', 0.0)
                    v1_lon2 = hdr.get('V1_LON', 0.0)
                    v2_lon2 = hdr.get('V2_LON', 0.0)
                    v3_lon2 = hdr.get('V3_LON', 0.0)

        # Calculate the value of SPICE parameters for the added files
        lon1 = np.array([v0_lon1, v1_lon1, v2_lon1, v3_lon1])
        lon2 = np.array([v0_lon2, v1_lon2, v2_lon2, v3_lon2])
        lat1 = np.array([v0_lat1, v1_lat1, v2_lat1, v3_lat1])
        lat2 = np.array([v0_lat2, v1_lat2, v2_lat2, v3_lat2])

        dy1 = lat1
        dx1 = lon1
        dy2 = lat2
        dx2 = lon2

        cx = np.zeros(4)
        cy = np.zeros(4)

        if max(dy1) > max(dy2):
            # First footprint is at top
            cx = [dx1[0], dx2[1], dx2[2], dx1[3]]
            cy = [dy1[0], dy2[1], dy2[2], dy1[3]]
        elif max(dy1) < max(dy2):
            # Second footprint is at top
            cx = [dx2[0], dx1[1], dx1[2], dx2[3]]
            cy = [dy2[0], dy1[1], dy1[2], dy2[3]]

        # Calculate mean values
        sat_alt_mean = np.mean(sat_alt_all)
        sol_ang_mean = np.mean(sol_ang_all)
        phs_ang_mean = np.mean(phs_ang_all)
        emi_ang_mean = np.mean(emi_ang_all)


        temperature_mean = np.mean(temperature_all)
        gain_mean = np.mean(gain_all)
        scd_filter_mean = np.mean(scd_filter_all)
        sat_lat_mean = np.mean(sat_lat_all)
        sat_lon_mean = np.mean(sat_lon_all)
        bore_lat_mean = np.mean(bore_lat_all)
        bore_lon_mean = np.mean(bore_lon_all)

        # Calculate average Local Solar Time
        lst_total_seconds = [hr * 3600 + mn * 60 + sc for hr, mn, sc in zip(lst_hr_all, lst_min_all, lst_sec_all)]
        lst_avg_seconds = np.mean(lst_total_seconds)
        lst_avg_seconds = lst_avg_seconds % 86400  # Ensure it's within 0 to 86400 seconds

        lst_avg_hr = int(lst_avg_seconds // 3600)
        lst_avg_min = int((lst_avg_seconds % 3600) // 60)
        lst_avg_sec = lst_avg_seconds % 60



        # Put all the SPICE values in a dictionary
        SPICE_val = {
            'sat_alt_mean': sat_alt_mean,
            'sol_ang_mean': sol_ang_mean,
            'phs_ang_mean': phs_ang_mean,
            'emi_ang_mean': emi_ang_mean,


            'temperature_mean': temperature_mean,
            'gain_mean': gain_mean,
            'scd_filter_mean': scd_filter_mean,
            'sat_lat_mean': sat_lat_mean,
            'sat_lon_mean': sat_lon_mean,
            'bore_lat_mean': bore_lat_mean,
            'bore_lon_mean': bore_lon_mean,

            'LST_HR': lst_avg_hr,
            'LST_MIN': lst_avg_min,
            'LST_SEC': lst_avg_sec,

            'v0_lat': cy[0],
            'v1_lat': cy[1],
            'v2_lat': cy[2],
            'v3_lat': cy[3],
            'v0_lon': cx[0],
            'v1_lon': cx[1],
            'v2_lon': cx[2],
            'v3_lon': cx[3]
        }

        print('')
        print('Total no. of files matching the start and end times provided by user:', count_ind_UTC_match)
        print('List of files added:')
        for f in files_to_add_arr:
            print(f)

        print('')
        print('Start UTC provided by user:', UTC_start_user)
        print('End UTC provided by user:', UTC_end_user)

        # Now write as a FITS file
        UTC_start_str = utc_to_timestr(UTC_start_actual)
        UTC_end_str = utc_to_timestr(UTC_end_actual)

        summed_SCD_opfname = os.path.join(op_subdir_path, f'ch2_cla_l1_{UTC_start_str}_{UTC_end_str}{fname_ext}')
        # Calculate the exposure time based on the number of files added - each file is 8s
        expotime = count_ind_UTC_match * 8.0
        print('Exposure time (seconds):', expotime)

        # Write a fits file
        write_summed_fits(summed_SCD_opfname, data_SCD_summed, ip_file, expotime, [UTC_start_actual], [UTC_end_actual], SPICE_val)
        print('')
        print('***************')
        print('Program CLASS_add_L1_files_time.py completed.')
        print('Output file name:', summed_SCD_opfname)
        print('Output files of this program are stored in the following location:')
        print(op_subdir_path)

        return(summed_SCD_opfname)

    else:
        print("Could not find files to add for the user specified start and end times.")
        print("Please try again with proper values.")

if __name__ == '__main__':
    # Example usage
    ipfile_dir = '/home/sammy/Downloads/Processing/M/20210828/flare_68/'
    UTC_start_user = '2021-08-28T06:22:57.451'
    UTC_end_user = '2021-08-28T06:23:29.451'
    op_dir_location = '/home/sammy/Downloads/Processing/M/20210828/flare_68/integrated/'
    class_add_l1_files_time(ipfile_dir, UTC_start_user, UTC_end_user, op_dir_location)
