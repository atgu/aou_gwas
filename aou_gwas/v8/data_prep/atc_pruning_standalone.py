#!/usr/bin/env python3
"""
ATC Code Pruning Script for Hail

This script analyzes ATC codes and their case/control overlap.
For ATC codes with >90% case overlap with their parent code, we keep only the parent.

This script directly processes data without function calls for easier understanding.
"""

import os
import pandas as pd
import numpy as np
import hail as hl
from google.cloud import storage
import logging
import time

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

#############################
# CONFIGURATION PARAMETERS
#############################

# Input and output paths
case_counts_file = "gs://aou_analysis/v8/data/phenotype/summary/r_drug_summary.ht"  # Path to Hail Table with case counts
output_path = "gs://aou_analysis/v8/data/phenotype/summary/"
overlap_threshold = 90  # percentage threshold for pruning

# Ensure the output path ends with /
if not output_path.endswith('/'):
    output_path += '/'

#############################
# INITIALIZE HAIL
#############################

logger.info(f"Initializing Hail")

# Initialize Hail with increased memory settings
hl.init(default_reference="GRCh38", 
        master='local[4]',
        app_name='aou_v8_atc_curation',
        driver_memory='highmem')

#############################
# READ CASE COUNTS FROM HAIL TABLE
#############################

logger.info(f"Reading ATC code case counts from Hail Table: {case_counts_file}")

# Read the Hail Table
counts_ht = hl.read_table(case_counts_file)
logger.info(f"Hail Table schema: {counts_ht.describe()}")

# Convert to pandas for easier processing
case_counts_df = counts_ht.to_pandas()
logger.info(f"Case counts dataframe has {len(case_counts_df)} rows")
logger.info(f"Case counts dataframe columns: {case_counts_df.columns}")

# Determine the correct column names
code_col = 'atc_code' if 'atc_code' in case_counts_df.columns else 'code'
count_col = 'case_count' if 'case_count' in case_counts_df.columns else 'n_cases'
logger.info(f"Using columns: {code_col} and {count_col}")

# Get all ATC codes and create a dictionary of case counts
atc_codes = case_counts_df[code_col].tolist()
case_counts = dict(zip(case_counts_df[code_col], case_counts_df[count_col]))

# Get total number of samples
total_samples_col = 'total_samples' if 'total_samples' in case_counts_df.columns else None
if total_samples_col:
    total_samples = case_counts_df[total_samples_col].iloc[0]
    logger.info(f"Total number of samples: {total_samples}")
else:
    # Estimate total samples as maximum case count + estimated controls
    estimated_max_cases = max(case_counts.values())
    total_samples = estimated_max_cases * 10  # Assuming roughly 10x more controls than max cases
    logger.info(f"Total number of samples estimated as: {total_samples}")

#############################
# CREATE RESULTS DATAFRAME
#############################

# Create result dataframe to store our analysis
results = pd.DataFrame({
    'atc_code': atc_codes,
    'case_count': [case_counts.get(code, 0) for code in atc_codes],
    'is_kept': True,  # Initially assume all codes are kept
    'closest_kept_parent': None,  # Will show the closest parent that is kept
    'pruned_by': None,  # Will store which ancestor code caused pruning
    'overlap_percentage': 0.0,  # Will store overlap with closest non-pruned ancestor
    'data_inconsistency': False  # Will be True if child has more cases than parent
})

#############################
# HIERARCHICAL PRUNING ALGORITHM
#############################

logger.info("Starting hierarchical pruning algorithm")

# First, sort codes by length so we process parents before children
results['code_length'] = results['atc_code'].astype(str).str.len()
results = results.sort_values(['code_length', 'atc_code']).reset_index(drop=True)

# Function to check if code2 is an ancestor of code1
def is_ancestor(code1, code2):
    return code1.startswith(code2) and code1 != code2

# Process codes in hierarchical order (shorter codes first)
for i, row in results.iterrows():
    current_code = row['atc_code']
    current_case_count = row['case_count']
    
    # Skip codes with no cases
    if current_case_count == 0:
        continue
    
    # Find closest non-pruned ancestor (if any)
    ancestor_codes = []
    for j, potential_ancestor in results.iloc[:i].iterrows():
        if potential_ancestor['is_kept'] and is_ancestor(current_code, potential_ancestor['atc_code']):
            ancestor_codes.append(potential_ancestor['atc_code'])
    
    # Sort ancestors by length (descending) to get the closest one
    ancestor_codes.sort(key=len, reverse=True)
    
    if ancestor_codes:
        closest_ancestor = ancestor_codes[0]
        results.loc[i, 'closest_kept_parent'] = closest_ancestor
        ancestor_case_count = case_counts.get(closest_ancestor, 0)
        
        # Data validation - a child should never have more cases than its parent
        if current_case_count > ancestor_case_count:
            logger.warning(f"Data inconsistency: {current_code} has {current_case_count} cases, "
                          f"which is more than its parent {closest_ancestor} with {ancestor_case_count} cases")
            
            # Flag this inconsistency in the results
            results.loc[i, 'data_inconsistency'] = True
            
            # In a proper hierarchy, we should never see this, but let's handle it gracefully
            # Since the parent should contain all children, we'll use the max value
            effective_ancestor_count = max(ancestor_case_count, current_case_count)
            
            if effective_ancestor_count == 0:
                overlap_percentage = 0
            else:
                # Calculate proper overlap - what percentage of the parent's cases does this child represent?
                overlap_percentage = (current_case_count / effective_ancestor_count) * 100
        else:
            # Normal case - Calculate overlap percentage - child cases as percentage of parent cases
            if ancestor_case_count > 0:
                overlap_percentage = (current_case_count / ancestor_case_count) * 100
            else:
                overlap_percentage = 0
        
        results.loc[i, 'overlap_percentage'] = overlap_percentage
        
        # Decide if we should prune - only prune if child contains >90% of parent's cases
        if overlap_percentage > overlap_threshold:
            results.loc[i, 'is_kept'] = False
            results.loc[i, 'pruned_by'] = closest_ancestor
            logger.info(f"Pruning {current_code} (overlap {overlap_percentage:.2f}% with ancestor {closest_ancestor})")
        else:
            logger.info(f"Keeping {current_code} (overlap {overlap_percentage:.2f}% with ancestor {closest_ancestor})")

# After initial pruning, make another pass to confirm all codes are correctly kept/pruned
corrected_count = 0
for i, row in results.iterrows():
    # Skip top-level codes and already pruned codes
    if row['closest_kept_parent'] is None or not row['is_kept']:
        continue
    
    current_code = row['atc_code']
    current_case_count = row['case_count']
    parent_code = row['closest_kept_parent']
    parent_case_count = case_counts.get(parent_code, 0)
    
    if parent_case_count > 0:
        overlap_percentage = (current_case_count / parent_case_count) * 100
        
        # If overlap is below threshold, the code should be kept
        if overlap_percentage <= overlap_threshold and not row['is_kept']:
            results.loc[i, 'is_kept'] = True
            results.loc[i, 'pruned_by'] = None
            logger.info(f"Corrected: Keeping {current_code} (overlap {overlap_percentage:.2f}% with parent {parent_code})")
            corrected_count += 1

if corrected_count > 0:
    logger.info(f"Corrected {corrected_count} codes that were incorrectly marked for pruning")

logger.info(f"Hierarchical pruning complete. Keeping {results['is_kept'].sum()} out of {len(results)} ATC codes")

#############################
# ANALYZE DIFFERENT CUTOFFS
#############################

# Check how many codes would be kept/removed at different thresholds
logger.info("Analyzing pruning results at different overlap percentage thresholds:")
thresholds = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]
threshold_results = []

# Create a copy of results to avoid modifying the original
threshold_analysis = results.copy()

for threshold in thresholds:
    # Reset is_kept for this threshold analysis
    threshold_analysis['is_kept'] = True
    
    # Re-run pruning with this threshold
    for i, row in threshold_analysis.iterrows():
        current_code = row['atc_code']
        
        # Skip codes with no cases or already marked as not kept
        if row['case_count'] == 0:
            continue
        
        # Find closest non-pruned ancestor at this threshold
        ancestor_codes = []
        for j, potential_ancestor in threshold_analysis.iloc[:i].iterrows():
            if potential_ancestor['is_kept'] and is_ancestor(current_code, potential_ancestor['atc_code']):
                ancestor_codes.append(potential_ancestor['atc_code'])
        
        # Sort ancestors by length (descending)
        ancestor_codes.sort(key=len, reverse=True)
        
        if ancestor_codes:
            closest_ancestor = ancestor_codes[0]
            ancestor_case_count = case_counts.get(closest_ancestor, 0)
            
            # Calculate overlap and decide whether to prune
            if ancestor_case_count > 0:
                overlap_percentage = (row['case_count'] / ancestor_case_count) * 100
                if overlap_percentage > threshold * 100:
                    threshold_analysis.loc[i, 'is_kept'] = False
    
    # Count kept codes at this threshold
    kept_at_threshold = threshold_analysis['is_kept'].sum()
    pruned_at_threshold = len(threshold_analysis) - kept_at_threshold
    
    # Log and store results
    logger.info(f"Threshold {threshold:.2f}: Keep {kept_at_threshold} codes, Remove {pruned_at_threshold} codes")
    threshold_results.append({
        'threshold': threshold,
        'kept': kept_at_threshold,
        'removed': pruned_at_threshold
    })

# Create DataFrame for threshold analysis results
threshold_df = pd.DataFrame(threshold_results)
logger.info(f"Threshold analysis complete. Using {overlap_threshold}% for final pruning.")

#############################
# PREPARE FINAL RESULTS
#############################

# Get codes to keep/remove based on the is_kept flag
codes_to_keep = results.loc[results['is_kept'], 'atc_code'].tolist()
codes_to_remove = results.loc[~results['is_kept'], 'atc_code'].tolist()

logger.info(f"After pruning: keeping {len(codes_to_keep)} out of {len(results)} ATC codes")

# Create final detailed results dataframe
final_results = results[['atc_code', 'case_count', 'is_kept', 'pruned_by', 'overlap_percentage', 
                         'closest_kept_parent', 'data_inconsistency']]
final_results.rename(columns={'is_kept': 'should_keep'}, inplace=True)

#############################
# COUNT FINAL CASES/CONTROLS
#############################

# Create DataFrame for final counts more efficiently
final_counts = pd.DataFrame({'atc_code': codes_to_keep})

# Use case_counts to get case counts for each kept code
final_counts['n_cases'] = final_counts['atc_code'].map(case_counts).fillna(0)

# Calculate controls as total - cases
final_counts['n_controls'] = total_samples - final_counts['n_cases']

#############################
# WRITE RESULTS TO GCS
#############################

logger.info("Writing results to GCS...")

# Helper function to write DataFrame to GCS
def write_df_to_gcs(df, gcs_path):
    # Parse the GCS path
    path_parts = gcs_path[5:].split('/', 1)
    bucket_name = path_parts[0]
    blob_name = path_parts[1]
    
    # Initialize GCS client
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    
    # Convert to CSV and upload
    csv_data = df.to_csv(index=False)
    blob = bucket.blob(blob_name if blob_name.endswith('.csv') else f"{blob_name}.csv")
    blob.upload_from_string(csv_data, content_type='text/csv')
    
    logger.info(f"CSV data written to: gs://{bucket_name}/{blob.name}")

# Write all results to GCS
write_df_to_gcs(final_results, f"{output_path}overlap_results")
write_df_to_gcs(pd.DataFrame({'atc_code': codes_to_keep}), f"{output_path}codes_to_keep")
write_df_to_gcs(final_counts, f"{output_path}final_counts")
write_df_to_gcs(threshold_df, f"{output_path}threshold_analysis")

logger.info("Analysis complete. Results written to GCS as CSV files.")
