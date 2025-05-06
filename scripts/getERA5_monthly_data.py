import cdsapi
import pandas as pd
import xarray as xr
from datetime import datetime
import os


import logging

logging.basicConfig(
    filename='/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/ERA5/DownloadMonthly.log',  # Ensure this matches your redirection target
    level=logging.DEBUG,       # Change based on your needs (DEBUG, INFO, WARNING, etc.)
    format='%(asctime)s - %(levelname)s - %(message)s'
)

logging.info("Script started.")

# Initialize the CDS API client
print("Initialize the CDS API client")
c = cdsapi.Client()

# Input and output file paths
input_csv = "/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/samples_europe_pass.csv"
output_csv = "/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/ERA5/ERA5_Europe_Monthly_Pass.csv"

# Read the input CSV file
print("Read the input CSV file")
df = pd.read_csv(input_csv)
try:
    df['date'] = pd.to_datetime(df['date'], format='%Y-%m-%d')
except KeyError:
    # If 'date' column is missing, extract date from 'sample' column and convert
    df['date'] = pd.to_datetime(df['sample'].str[-10:], format='%Y-%m-%d', errors='coerce')


# List of variables to retrieve for each dataset
#variables = [
#    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
#    '2m_dewpoint_temperature', 'mean_sea_level_pressure', 'surface_pressure',
#    'total_precipitation', 'surface_solar_radiation_downwards',
#    'surface_thermal_radiation_downwards', 'soil_temperature_level_1',
#    'soil_temperature_level_2', 'soil_temperature_level_3',
#    'soil_temperature_level_4', 'volumetric_soil_water_layer_1',
#    'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3',
#    'volumetric_soil_water_layer_4', 'skin_temperature', 'snow_depth',
#    'sea_surface_temperature', 'mean_wave_direction', 'mean_wave_period',
#    'significant_height_of_combined_wind_waves_and_swell'
#]

variables = ['10m_u_component_of_wind', '10m_v_component_of_wind',
             "2m_dewpoint_temperature",
            '2m_dewpoint_temperature', 'mean_sea_level_pressure',
            "2m_temperature", 
            "surface_pressure",
            "total_precipitation",
            'surface_solar_radiation_downwards',
            'surface_thermal_radiation_downwards',
            "maximum_2m_temperature_since_previous_post_processing",
            "minimum_2m_temperature_since_previous_post_processing",
            "mean_snowfall_rate",
            "precipitation_type",
            "snowfall",
            "soil_temperature_level_1",'soil_temperature_level_2', 'soil_temperature_level_3',
            'soil_temperature_level_4', 'volumetric_soil_water_layer_1',
            'volumetric_soil_water_layer_2',
            "soil_type",
            "high_vegetation_cover",
            "leaf_area_index_high_vegetation",
            "leaf_area_index_low_vegetation",
            "low_vegetation_cover", 'skin_temperature', 'snow_depth',
            'sea_surface_temperature', 'mean_wave_direction', 'mean_wave_period',
            'significant_height_of_combined_wind_waves_and_swell']


# Define available datasets
datasets = {
    'reanalysis-era5-single-levels-monthly-means': 'monthly'
}

#datasets = {
#    'reanalysis-era5-single-levels': 'daily'
#}
#

# Create columns in the DataFrame for each variable and dataset
print("# Create columns in the DataFrame for each variable and dataset")
for dataset, frequency in datasets.items():
    for variable in variables:
        df[f"{dataset}_{variable}_{frequency}"] = None

# Function to download and retrieve a specific variable from ERA5
#def get_era5_point_data(dataset, variable, year, month, day, latitude, longitude):
#    temp_file = f"temp_{dataset}_{variable}_{latitude}_{longitude}_{year}{month:02d}{day if day else ''}.nc"
#    print(temp_file)
#    c.retrieve(
#        dataset,
#        {
#            'product_type': 'reanalysis',
#            'format': 'netcdf',
#            'variable': variable,
#            'year': str(year),
#            'month': [f"{month:02d}"],
#            'day': [f"{day:02d}"] if day else None,
#            'time': '00:00',
#            'area': [latitude, longitude, latitude, longitude],
#        },
#        temp_file
#    )
#    # Open the NetCDF file to extract the variable value
#    ## THIS part needs to be rewritten  in order to get the variable names correctly and parse them correctly wiht ''
#    with xr.open_dataset(temp_file) as ds:
#        var = list(ds.data_vars.items())
#        varname, rest = var[0]
#        value = ds[varname].values[0, 0, 0]  # Assumes time, lat, lon dimensions
#    if os.path.exists(temp_file):
#        os.remove(temp_file)
#    return value

def get_era5_point_data(dataset, variable, year, month, day, latitude, longitude):
    temp_file = f"temp_{dataset}_{variable}_{latitude}_{longitude}_{year}{month:02d}{day if day else ''}.nc"
    print(f"Temporary file: {temp_file}")
    try:
        # Request data from the CDS API
        c.retrieve(
            dataset,
            {
                'product_type': ["monthly_averaged_reanalysis"],
                'data_format': 'netcdf',
                'variable': variable,
                'year': str(year),
                'month': [f"{month:02d}"],
                'day': [f"{day:02d}"] if day else None,
                'time': '00:00',
                'area': [latitude, longitude, latitude, longitude],
            },
            temp_file
        )      
        # Open the NetCDF file and extract the data variable value
        with xr.open_dataset(temp_file) as ds:
            # Dynamically fetch the first variable name
            varname = next(iter(ds.data_vars))
            print(f"Extracted variable: {varname}")
            # Retrieve the value assuming dimensions (time, latitude, longitude)
            value = ds[varname].values[0]
            print("AVERAGE VALUE:", value)
    except Exception as e:
        print(f"An error occurred: {e}")
        value = None  # Return None or handle differently in case of failure 
    finally:
        # Ensure the temporary file is deleted even if an error occurs
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                print(f"Temporary file {temp_file} deleted successfully.")
            except Exception as cleanup_error:
                print(f"Failed to delete temporary file {temp_file}: {cleanup_error}")
    return value

print("Starting Download")
# Iterate over each row and download data for all datasets and variables
for idx, row in df.iterrows():
    latitude, longitude, date, sample = row['latitude'], row['longitude'], row['date'], row['sample']
    year, month, day = date.year, date.month, date.day
    print(f"Collecting data for {sample}.")
    for dataset, frequency in datasets.items():
        for variable in variables:
            try:
                if frequency == 'daily':
                    value = get_era5_point_data(dataset, variable, year, month, day, latitude, longitude)
                    df.at[idx, f"{dataset}_{variable}_daily"] = value
                elif frequency == 'monthly':
                    value = get_era5_point_data(dataset, variable, year, month, None, latitude, longitude)
                    df.at[idx, f"{dataset}_{variable}_monthly"] = value
            except Exception as e:
                print(f"Failed to retrieve {variable} for {date} at {latitude}, {longitude}: {e}")

# Save the extended DataFrame to a CSV file
df.to_csv(output_csv, index=False)
print(f"Data with all ERA5 datasets saved to {output_csv}")
