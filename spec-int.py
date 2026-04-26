#!/usr/bin/env python3
"""
This script calculates the integral of spectral data between two wavelength points
with baseline correction by connecting the endpoints with a straight line.

The script is designed to work directly with raw output files from spectrometer instruments,
handling metadata rows and finding the actual spectral data automatically.

Usage:
  python spec-int.py [input_file] [start_wavelength] [end_wavelength] [column_names]

Arguments:
  input_file - Path to the CSV file (default: example.csv)
  start_wavelength - Starting wavelength for integration (default: 400)
  end_wavelength - Ending wavelength for integration (default: 650)
  column_names - Names of intensity columns to process, comma-separated (default: all columns)
                 Use "all" to process all non-empty data columns

Examples:
  python spec-int.py data.csv 370 650
  python spec-int.py data.csv 300 500 H1,H3,H5
  python spec-int.py data.csv 250 600 all
"""

import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Union, Optional, Set

class SpectrometerDataProcessor:
    """Class to handle spectrometer data files and perform integration."""
    
    def __init__(self, filename: str):
        """Initialize with a spectrometer data file."""
        self.filename = filename    
        self.wavelength_column = None
        self.headers = []
        self.data_columns = []
        self.metadata = {}
        self.wavelengths = []
        self.raw_data = []
        self.data_by_column = {}
        
    def read_file(self) -> None:
        """
        Read the spectrometer data file, handling metadata and identifying data structure.
        """
        print(f"Reading data from '{self.filename}'...")
        
        # Read all lines from the file
        with open(self.filename, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # Find the header row (typically contains "Wavelength")
        header_row_idx = -1
        for i, line in enumerate(lines):
            if "Wavelength" in line:
                header_row_idx = i
                break
        
        if header_row_idx == -1:
            raise ValueError("Could not find wavelength column in the file. Is this a valid spectrometer data file?")
        
        # Extract metadata from lines before the header
        for i in range(header_row_idx):
            if lines[i] and ',' in lines[i]:
                parts = lines[i].split(',', 1)
                if parts[0].strip():
                    self.metadata[parts[0].strip()] = parts[1].strip()
        
        # Parse the header row
        headers = lines[header_row_idx].split(',')
        
        # Find the wavelength column
        wavelength_idx = -1
        for i, header in enumerate(headers):
            if header.strip() == "Wavelength":
                wavelength_idx = i
                break
        
        if wavelength_idx == -1:
            raise ValueError("Could not find wavelength column in the header row.")
        
        self.wavelength_column = wavelength_idx
        self.headers = headers
        
        # Identify data columns (non-empty columns after wavelength)
        self.data_columns = []
        for i in range(wavelength_idx + 1, len(headers)):
            if headers[i].strip():
                self.data_columns.append((i, headers[i].strip()))
        
        if not self.data_columns:
            raise ValueError("No data columns found in the file.")
        
        # Read data rows (rows after the header that have a valid wavelength)
        data_rows = []
        for i in range(header_row_idx + 1, len(lines)):
            row = lines[i].split(',')
            if len(row) > wavelength_idx and row[wavelength_idx].strip().isdigit():
                data_rows.append(row)
        
        if not data_rows:
            raise ValueError("No valid data rows found in the file.")
        
        self.raw_data = data_rows
        self.wavelengths = [int(row[wavelength_idx]) for row in data_rows]
        
        print(f"Found {len(self.data_columns)} data columns: {', '.join(col[1] for col in self.data_columns)}")
        print(f"Found {len(self.wavelengths)} wavelength points from {min(self.wavelengths)} to {max(self.wavelengths)} nm")
        
    def extract_column_data(self, selected_columns: Optional[List[str]] = None) -> Dict[str, List[Dict]]:
        """
        Extract data for selected columns and organize it for processing.
        
        Args:
            selected_columns: List of column names to extract, or None for all columns
            
        Returns:
            Dictionary mapping column names to lists of wavelength/intensity data points
        """
        # Filter columns based on selection
        columns_to_process = []
        
        if not selected_columns or selected_columns == ["all"]:
            # Process all data columns
            columns_to_process = self.data_columns
        else:
            # Process only selected columns
            column_indices = {col_name: col_idx for col_idx, col_name in self.data_columns}
            for col_name in selected_columns:
                if col_name in column_indices:
                    columns_to_process.append((column_indices[col_name], col_name))
                else:
                    print(f"Warning: Column '{col_name}' not found in the data file.")
        
        if not columns_to_process:
            raise ValueError("No valid columns to process.")
        
        # Extract data by column
        data_by_column = {}
        
        for col_idx, col_name in columns_to_process:
            data_points = []
            
            for row_idx, row in enumerate(self.raw_data):
                if len(row) <= col_idx:
                    # Skip if row doesn't have enough columns
                    continue
                
                wavelength = self.wavelengths[row_idx]
                value_str = row[col_idx].strip()
                
                # Handle non-numeric values
                try:
                    intensity = float(value_str) if value_str else None
                except ValueError:
                    intensity = value_str  # Keep as string if it's not a number
                
                data_points.append({
                    'Wavelength': wavelength,
                    'Intensity': intensity
                })
            
            # Sort by wavelength just to be safe
            data_points.sort(key=lambda x: x['Wavelength'])
            
            # Filter out None values
            data_points = [p for p in data_points if p['Intensity'] is not None]
            
            if data_points:
                data_by_column[col_name] = data_points
        
        return data_by_column

def find_nearest_wavelength_indices(data: List[Dict], target_wavelengths: List[int]) -> List[int]:
    """
    Find the indices of the wavelengths closest to the target wavelengths.
    
    Args:
        data: List of dictionaries containing wavelength and intensity data
        target_wavelengths: List of wavelengths to find
        
    Returns:
        List of indices in the data that correspond to the nearest wavelengths
    """
    indices = []
    wavelengths = [point['Wavelength'] for point in data]
    
    for target in target_wavelengths:
        nearest_idx = min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i] - target))
        indices.append(nearest_idx)
        
    return indices

def calculate_baseline_corrected_integral(data: List[Dict], start_idx: int, end_idx: int) -> Tuple[float, Dict]:
    """
    Calculate the integral between two wavelengths with baseline correction.
    
    Args:
        data: List of dictionaries containing wavelength and intensity data
        start_idx: Index of the starting wavelength
        end_idx: Index of the ending wavelength
        
    Returns:
        Tuple containing the calculated integral value and info about any issues
    """
    start_point = data[start_idx]
    end_point = data[end_idx]
    
    info = {
        'overflow_points': 0,
        'start_wavelength': start_point['Wavelength'],
        'end_wavelength': end_point['Wavelength'],
        'warnings': []
    }
    
    # Check if endpoints contain non-numeric values
    if isinstance(start_point['Intensity'], str):
        info['warnings'].append(f"Start point at wavelength {start_point['Wavelength']} has overflow value '{start_point['Intensity']}'. Using nearby valid data.")
        # Find the nearest valid data point after the start
        for i in range(start_idx + 1, end_idx):
            if isinstance(data[i]['Intensity'], (int, float)):
                start_idx = i
                start_point = data[start_idx]
                info['adjusted_start_wavelength'] = start_point['Wavelength']
                break
    
    if isinstance(end_point['Intensity'], str):
        info['warnings'].append(f"End point at wavelength {end_point['Wavelength']} has overflow value '{end_point['Intensity']}'. Using nearby valid data.")
        # Find the nearest valid data point before the end
        for i in range(end_idx - 1, start_idx, -1):
            if isinstance(data[i]['Intensity'], (int, float)):
                end_idx = i
                end_point = data[end_idx]
                info['adjusted_end_wavelength'] = end_point['Wavelength']
                break
    
    # Check if we have valid numeric values at endpoints
    if isinstance(start_point['Intensity'], str) or isinstance(end_point['Intensity'], str):
        info['warnings'].append("Could not find valid numeric values at endpoints. Cannot calculate integral.")
        return 0.0, info
    
    # Calculate the baseline (straight line connecting the endpoints)
    x1, y1 = start_point['Wavelength'], start_point['Intensity']
    x2, y2 = end_point['Wavelength'], end_point['Intensity']
    
    # Straight line equation: y = mx + b
    m = (y2 - y1) / (x2 - x1) if x2 != x1 else 0
    b = y1 - m * x1
    
    # Calculate the integral
    integral = 0.0
    for i in range(start_idx, end_idx):
        current = data[i]
        next_point = data[i + 1]
        
        # Skip if either point is not numeric
        if (not isinstance(current['Intensity'], (int, float)) or 
            not isinstance(next_point['Intensity'], (int, float))):
            info['overflow_points'] += 1
            continue
        
        x_current, y_current = current['Wavelength'], current['Intensity']
        x_next, y_next = next_point['Wavelength'], next_point['Intensity']
        
        # Calculate baseline values at these wavelengths
        baseline_current = m * x_current + b
        baseline_next = m * x_next + b
        
        # Calculate baseline-corrected areas using trapezoidal rule
        corrected_current = y_current - baseline_current
        corrected_next = y_next - baseline_next
        
        # Area of the trapezoid
        area = 0.5 * (corrected_current + corrected_next) * (x_next - x_current)
        integral += area
    
    return integral, info

def plot_spectrum_with_integration(data: List[Dict], start_idx: int, end_idx: int, column_name: str, integral: float) -> str:
    """
    Create a plot of the spectrum showing the integration area and baseline correction.
    
    Args:
        data: List of dictionaries containing wavelength and intensity data
        start_idx: Index of the starting wavelength
        end_idx: Index of the ending wavelength
        column_name: Name of the column being plotted
        integral: Calculated integral value
        
    Returns:
        Filename of the saved plot
    """
    wavelengths = [point['Wavelength'] for point in data]
    intensities = [point['Intensity'] if isinstance(point['Intensity'], (int, float)) else np.nan for point in data]
    
    # Convert to numpy arrays for easier manipulation
    wavelengths = np.array(wavelengths)
    intensities = np.array(intensities)
    
    start_wavelength = data[start_idx]['Wavelength']
    end_wavelength = data[end_idx]['Wavelength']
    start_intensity = data[start_idx]['Intensity']
    end_intensity = data[end_idx]['Intensity']
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Plot the full spectrum
    plt.plot(wavelengths, intensities, 'b-', label='Spectrum')
    
    # Highlight the integration region
    plt.fill_between(
        wavelengths[start_idx:end_idx+1],
        intensities[start_idx:end_idx+1],
        np.interp(wavelengths[start_idx:end_idx+1], [start_wavelength, end_wavelength], [start_intensity, end_intensity]),
        alpha=0.3,
        color='green',
        label='Integrated Area'
    )
    
    # Draw the baseline
    plt.plot([start_wavelength, end_wavelength], [start_intensity, end_intensity], 'r--', label='Baseline')
    
    # Add markers for the integration limits
    plt.plot([start_wavelength], [start_intensity], 'ro', markersize=8)
    plt.plot([end_wavelength], [end_intensity], 'ro', markersize=8)
    
    # Annotate the plot
    plt.title(f'Baseline-corrected Integration for {column_name}\nIntegral = {integral:.4f}')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save the plot
    output_filename = f"integration_{column_name}_{start_wavelength}-{end_wavelength}.png"
    plt.savefig(output_filename, dpi=150)
    plt.close()
    
    return output_filename

def main():
    """Main function to process command line arguments and calculate the integral."""
    # Set default values
    input_file = 'example.csv'
    start_wavelength = 400
    end_wavelength = 650
    selected_columns = ["all"]  # Process all columns by default
    
    # Process command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        try:
            start_wavelength = int(sys.argv[2])
        except ValueError:
            print(f"Error: Start wavelength '{sys.argv[2]}' is not a valid number. Using default: {start_wavelength}")
    if len(sys.argv) > 3:
        try:
            end_wavelength = int(sys.argv[3])
        except ValueError:
            print(f"Error: End wavelength '{sys.argv[3]}' is not a valid number. Using default: {end_wavelength}")
    if len(sys.argv) > 4:
        selected_columns = sys.argv[4].split(',')
    
    # Ensure start_wavelength is less than end_wavelength
    if start_wavelength > end_wavelength:
        start_wavelength, end_wavelength = end_wavelength, start_wavelength
        print(f"Warning: Start wavelength was greater than end wavelength. Swapped to {start_wavelength}-{end_wavelength}.")
    
    try:
        # Check if matplotlib is available for plotting
        has_matplotlib = True
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            has_matplotlib = False
            print("Warning: matplotlib not available. Plots will not be generated.")
        
        # Create processor and read data
        processor = SpectrometerDataProcessor(input_file)
        processor.read_file()
        
        # Print some metadata if available
        if 'Date' in processor.metadata:
            print(f"Measurement date: {processor.metadata['Date']}")
        if 'Time' in processor.metadata:
            print(f"Measurement time: {processor.metadata['Time']}")
        if 'Reader Type:' in processor.metadata:
            print(f"Instrument: {processor.metadata['Reader Type:']}")
        
        # Extract data for the selected columns
        data_by_column = processor.extract_column_data(selected_columns)
        
        processed_columns = list(data_by_column.keys())
        if not processed_columns:
            print("No columns were processed. Please check your column selection.")
            sys.exit(1)
            
        print(f"Processing columns: {', '.join(processed_columns)}")
        
        # Store results for all columns
        results = {}
        plot_files = []
        
        # Process each column
        for column in processed_columns:
            data = data_by_column[column]
            
            # Find the indices closest to the target wavelengths
            start_idx, end_idx = find_nearest_wavelength_indices(data, [start_wavelength, end_wavelength])
            
            # Calculate the baseline-corrected integral
            integral, info = calculate_baseline_corrected_integral(data, start_idx, end_idx)
            
            # Store results
            results[column] = {
                'integral': integral,
                'info': info
            }
            
            # Create a plot if matplotlib is available
            if has_matplotlib:
                try:
                    plot_file = plot_spectrum_with_integration(data, start_idx, end_idx, column, integral)
                    plot_files.append(plot_file)
                    print(f"Created plot: {plot_file}")
                except Exception as e:
                    print(f"Warning: Could not create plot for {column}: {str(e)}")
        
        # Print the results in a table format
        print("\nIntegration Results:")
        print("-" * 80)
        print(f"{'Column':<10} | {'Range':<20} | {'Integral':<15} | {'Warnings':<30}")
        print("-" * 80)
        
        for column in processed_columns:
            result = results[column]
            integral = result['integral']
            info = result['info']
            
            range_str = f"{info['start_wavelength']}-{info['end_wavelength']} nm"
            if 'adjusted_start_wavelength' in info or 'adjusted_end_wavelength' in info:
                range_str += " (adjusted)"
                
            warnings_str = ""
            if info['overflow_points'] > 0:
                warnings_str += f"{info['overflow_points']} overflow points"
                
            print(f"{column:<10} | {range_str:<20} | {integral:<15.4f} | {warnings_str:<30}")
            
        print("-" * 80)
            
    except FileNotFoundError:
        print(f"Error: Could not find file '{input_file}'")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()