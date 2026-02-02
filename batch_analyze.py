#!/usr/bin/env python3
"""
Batch Crown Analysis Script
===========================

Process multiple crown pairs and export results to a single CSV file.

Usage:
    python batch_analyze.py --input data_folder --output results.csv

The input folder should contain pairs of STL files with naming convention:
    - {specimen_id}_baseline.stl
    - {specimen_id}_adjusted.stl

Or use a CSV manifest file specifying the pairs:
    python batch_analyze.py --manifest pairs.csv --output results.csv

Manifest CSV format:
    specimen_id,baseline_path,adjusted_path,group
    1,baseline_1.stl,adjusted_1.stl,ExoCAD_AI
    2,baseline_2.stl,adjusted_2.stl,Human_Expert
"""

import argparse
import os
import sys
import csv
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd

# Import the analyzer from main module
from crown_analyzer import MeshAnalyzer


def find_stl_pairs(folder: str) -> List[Tuple[str, str, str]]:
    """
    Find matching baseline/adjusted STL pairs in a folder.
    
    Expected naming convention:
    - {id}_baseline.stl / {id}_adjusted.stl
    - {id}_pre.stl / {id}_post.stl
    - {id}_before.stl / {id}_after.stl
    
    Returns:
        List of tuples: (specimen_id, baseline_path, adjusted_path)
    """
    folder_path = Path(folder)
    stl_files = list(folder_path.glob("*.stl")) + list(folder_path.glob("*.STL"))
    
    # Group files by potential specimen ID
    pairs = []
    
    # Try different naming conventions
    baseline_suffixes = ['_baseline', '_pre', '_before', '_original', '_mill', '_milled']
    adjusted_suffixes = ['_adjusted', '_post', '_after', '_modified', '_adj']
    
    processed = set()
    
    for stl_file in stl_files:
        if stl_file in processed:
            continue
            
        stem = stl_file.stem.lower()
        
        # Check if this is a baseline file
        for b_suffix in baseline_suffixes:
            if stem.endswith(b_suffix):
                specimen_id = stem[:-len(b_suffix)]
                
                # Look for matching adjusted file
                for a_suffix in adjusted_suffixes:
                    for ext in ['.stl', '.STL']:
                        adjusted_name = specimen_id + a_suffix + ext
                        adjusted_path = folder_path / adjusted_name
                        if adjusted_path.exists():
                            pairs.append((
                                specimen_id,
                                str(stl_file),
                                str(adjusted_path)
                            ))
                            processed.add(stl_file)
                            processed.add(adjusted_path)
                            break
                    if stl_file in processed:
                        break
                break
    
    return pairs


def load_manifest(manifest_path: str) -> List[Dict]:
    """
    Load specimen pairs from a CSV manifest file.
    
    Expected columns: specimen_id, baseline_path, adjusted_path, [group]
    """
    pairs = []
    
    with open(manifest_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pairs.append({
                'specimen_id': row.get('specimen_id', row.get('id', '')),
                'baseline_path': row['baseline_path'],
                'adjusted_path': row['adjusted_path'],
                'group': row.get('group', 'Unknown')
            })
    
    return pairs


def batch_analyze(pairs: List[Dict], output_path: str, verbose: bool = True) -> pd.DataFrame:
    """
    Analyze multiple crown pairs and save results.
    
    Parameters
    ----------
    pairs : List[Dict]
        List of dictionaries with keys: specimen_id, baseline_path, adjusted_path, [group]
    output_path : str
        Path for output CSV file
    verbose : bool
        Print progress information
        
    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with all analysis results
    """
    analyzer = MeshAnalyzer()
    results_list = []
    
    total = len(pairs)
    
    for i, pair in enumerate(pairs):
        specimen_id = pair.get('specimen_id', f'specimen_{i+1}')
        group = pair.get('group', 'Unknown')
        baseline_path = pair['baseline_path']
        adjusted_path = pair['adjusted_path']
        
        if verbose:
            print(f"\n[{i+1}/{total}] Analyzing: {specimen_id} ({group})")
            print(f"  Baseline: {os.path.basename(baseline_path)}")
            print(f"  Adjusted: {os.path.basename(adjusted_path)}")
        
        try:
            # Verify files exist
            if not os.path.exists(baseline_path):
                print(f"  ERROR: Baseline file not found: {baseline_path}")
                continue
            if not os.path.exists(adjusted_path):
                print(f"  ERROR: Adjusted file not found: {adjusted_path}")
                continue
            
            # Run analysis
            results = analyzer.analyze(baseline_path, adjusted_path)
            
            vol = results['volume_analysis']
            dev = results['deviation_analysis']
            
            # Compile results row
            row = {
                'specimen_id': specimen_id,
                'group': group,
                'baseline_file': os.path.basename(baseline_path),
                'adjusted_file': os.path.basename(adjusted_path),
                'baseline_volume_mm3': vol['baseline_volume_mm3'],
                'adjusted_volume_mm3': vol['adjusted_volume_mm3'],
                'delta_volume_mm3': vol['delta_volume_mm3'],
                'percent_volume_removed': vol['percent_volume_removed'],
                'rms_um': dev['rms_um'],
                'mean_deviation_um': dev['mean_um'],
                'std_deviation_um': dev['std_um'],
                'median_deviation_um': dev['median_um'],
                'max_deviation_um': dev['max_deviation_um'],
                'min_deviation_um': dev['min_deviation_um'],
                'adjustment_area_mm2': dev['adjustment_area_mm2'],
                'percent_surface_adjusted': dev['percent_surface_adjusted'],
                'n_vertices_baseline': results['baseline']['vertices'],
                'n_vertices_adjusted': results['adjusted']['vertices'],
                'baseline_watertight': results['baseline']['is_watertight'],
                'adjusted_watertight': results['adjusted']['is_watertight'],
                'analysis_timestamp': results['timestamp'],
                'status': 'Success'
            }
            
            results_list.append(row)
            
            if verbose:
                print(f"  ✓ Delta Volume: {vol['delta_volume_mm3']:.4f} mm³")
                print(f"  ✓ RMS Deviation: {dev['rms_um']:.2f} µm")
                
        except Exception as e:
            print(f"  ERROR: Analysis failed - {str(e)}")
            results_list.append({
                'specimen_id': specimen_id,
                'group': group,
                'baseline_file': os.path.basename(baseline_path),
                'adjusted_file': os.path.basename(adjusted_path),
                'status': f'Error: {str(e)}'
            })
    
    # Create DataFrame
    results_df = pd.DataFrame(results_list)
    
    # Save to CSV
    results_df.to_csv(output_path, index=False)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Batch analysis complete!")
        print(f"Results saved to: {output_path}")
        print(f"Total specimens: {total}")
        print(f"Successful: {len([r for r in results_list if r.get('status') == 'Success'])}")
        print(f"Failed: {len([r for r in results_list if r.get('status', '').startswith('Error')])}")
        
        # Print summary statistics by group
        if 'group' in results_df.columns and 'delta_volume_mm3' in results_df.columns:
            success_df = results_df[results_df['status'] == 'Success']
            if not success_df.empty:
                print(f"\n{'='*60}")
                print("Summary by Group:")
                print("-"*60)
                summary = success_df.groupby('group').agg({
                    'delta_volume_mm3': ['mean', 'std', 'count'],
                    'rms_um': ['mean', 'std']
                }).round(4)
                print(summary.to_string())
    
    return results_df


def create_sample_manifest(output_path: str):
    """Create a sample manifest CSV file."""
    sample_data = """specimen_id,baseline_path,adjusted_path,group
1,data/specimen_1_baseline.stl,data/specimen_1_adjusted.stl,ExoCAD_AI
2,data/specimen_2_baseline.stl,data/specimen_2_adjusted.stl,ExoCAD_AI
3,data/specimen_3_baseline.stl,data/specimen_3_adjusted.stl,Human_Expert
4,data/specimen_4_baseline.stl,data/specimen_4_adjusted.stl,Human_Expert
5,data/specimen_5_baseline.stl,data/specimen_5_adjusted.stl,3Shape_Automate
"""
    with open(output_path, 'w') as f:
        f.write(sample_data)
    print(f"Sample manifest created: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Batch analyze crown adjustment STL pairs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--input', '-i',
        help='Input folder containing STL pairs'
    )
    parser.add_argument(
        '--manifest', '-m',
        help='CSV manifest file specifying STL pairs'
    )
    parser.add_argument(
        '--output', '-o',
        default=f'batch_results_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv',
        help='Output CSV file path'
    )
    parser.add_argument(
        '--create-manifest',
        metavar='PATH',
        help='Create a sample manifest CSV file'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress output'
    )
    
    args = parser.parse_args()
    
    # Create sample manifest if requested
    if args.create_manifest:
        create_sample_manifest(args.create_manifest)
        return
    
    # Validate inputs
    if not args.input and not args.manifest:
        parser.print_help()
        print("\nError: Please specify either --input folder or --manifest file")
        sys.exit(1)
    
    # Build pairs list
    pairs = []
    
    if args.manifest:
        if not os.path.exists(args.manifest):
            print(f"Error: Manifest file not found: {args.manifest}")
            sys.exit(1)
        pairs = load_manifest(args.manifest)
        
    elif args.input:
        if not os.path.isdir(args.input):
            print(f"Error: Input folder not found: {args.input}")
            sys.exit(1)
        found_pairs = find_stl_pairs(args.input)
        pairs = [
            {'specimen_id': p[0], 'baseline_path': p[1], 'adjusted_path': p[2], 'group': 'Unknown'}
            for p in found_pairs
        ]
    
    if not pairs:
        print("Error: No STL pairs found!")
        sys.exit(1)
    
    print(f"Found {len(pairs)} STL pairs to analyze")
    
    # Run batch analysis
    batch_analyze(pairs, args.output, verbose=not args.quiet)


if __name__ == "__main__":
    main()
