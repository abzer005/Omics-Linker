import streamlit as st
from .common import *  # Importing common functionalities from the 'common' module
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
#import cupy as cp

def combine_dataframes(df1, df2):

    # Align columns by names
    common_columns = df1.columns.intersection(df2.columns)

    # Subset DataFrames to common columns
    df1_aligned = df1[common_columns]
    df2_aligned = df2[common_columns]

    # Combine the DataFrames
    combined = pd.concat([df1_aligned, df2_aligned], axis=0)

    # Reset index and add 'ID' column
    combined.reset_index(inplace=True)
    combined.rename(columns={"index": "ID"}, inplace=True)
    combined.set_index("ID", inplace=True)

    return combined


# def calculate_correlations_gpu(df, metabolome_ft, genome_ft):
#     """
#     Faster correlation calculation using GPUs.
#     """
#     transposed_df = df.T
#     length_metabolome = metabolome_ft.shape[0]
#     length_genome = genome_ft.shape[0]

#     # Convert data to CuPy arrays
#     metabolomics = cp.asarray(transposed_df.iloc[:, :length_metabolome].values)
#     asvs = cp.asarray(transposed_df.iloc[:, length_metabolome:length_metabolome + length_genome].values)

#     # Standardize data
#     metabolomics = (metabolomics - cp.mean(metabolomics, axis=0)) / cp.std(metabolomics, axis=0)
#     asvs = (asvs - cp.mean(asvs, axis=0)) / cp.std(asvs, axis=0)

#     # Calculate correlation matrix
#     correlation_matrix = cp.dot(asvs.T, metabolomics) / (asvs.shape[0] - 1)
#     return cp.asnumpy(correlation_matrix)  # Convert back to NumPy


# def calculate_metabolite_asv_correlations(df, metabolome_ft, genome_ft):

#     transposed_df = df.T
#     length_metabolome = metabolome_ft.shape[0]
#     length_genome = genome_ft.shape[0]

#     # Initialize an empty list to store results
#     correlation_results = []

#     # Define ranges for metabolites and ASVs
#     metabolite_columns = range(length_metabolome)
#     asv_columns = range(length_genome, length_metabolome + length_genome)

#     # Loop over each ASV column
#     for asv_index in asv_columns:
#         asv_column = transposed_df.iloc[:, asv_index]

#         # Calculate correlations for the current ASV column with all metabolite columns
#         correlations = np.array([
#             pearsonr(asv_column, transposed_df.iloc[:, metabolite_index])
#             for metabolite_index in metabolite_columns
#         ])

#         # Extract estimates and p-values into a matrix
#         estimates = correlations[:, 0]  # Correlation coefficients
#         p_values = correlations[:, 1]  # P-values
#         result_matrix = np.column_stack((estimates, p_values))

#         # Append the matrix to the list
#         correlation_results.append(result_matrix)


#     return transposed_df, correlation_results


def calculate_single_asv(asv_index, metabolomics, asvs):
    """
    Compute correlations for a single ASV column.
    """
    correlations = np.array([pearsonr(asvs[:, asv_index], metabolomics[:, j]) 
                     for j in range(metabolomics.shape[1])
                     ])
    
    # Extract correlation coefficients and p-values
    correlation_coefficients = correlations[:, 0]
    p_values = correlations[:, 1]

    # Calculate R^2 values
    r_squared = correlation_coefficients ** 2
    
    # Apply FDR correction to p-values
    _, fdr_corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    # Combine results into a single array
    result = np.column_stack((correlation_coefficients, p_values, fdr_corrected_p_values, r_squared))

    return result
    

def calculate_correlations_parallel(df, metabolome_ft, genome_ft):
    """
    Faster correlation calculation using parallel processing.
    """

    transposed_df = df.T
    length_metabolome = metabolome_ft.shape[0]
    length_genome = genome_ft.shape[0]

    metabolome_names = metabolome_ft.index
    asv_names = genome_ft.index

    # Extract metabolomics and ASV columns
    metabolomics = transposed_df.iloc[:, :length_metabolome].values
    asvs = transposed_df.iloc[:, length_metabolome:length_metabolome + length_genome].values

    # Use multiprocessing to calculate correlations in parallel
    with Pool() as pool:
        results = pool.starmap(calculate_single_asv, 
                               [(i, metabolomics, asvs) for i in range(asvs.shape[1])])
        
    results_with_indices = {
        asv_names[i]: pd.DataFrame(
            results[i],
            index=metabolome_names,
            columns=["Estimate", "P-value", "BH-Corrected P-Value", "R2"]
        )
        for i in range(len(asv_names))
    }

    return results_with_indices

def merge_asv_correlation_results(results, target_df, genome_ft):
    """
    Merge ASV correlation results into a dictionary with dataframes having the same rows as target_df.

    Parameters:
    results (dict): Dictionary of correlation results for each ASV.
    target_df (pd.DataFrame): The combined metabolomics and genomics dataframe.
    gen_ft (pd.DataFrame): Genomics feature table.

    Returns:
    dict: A dictionary of merged dataframes for each ASV.
    """
    merged_list = {}

    empty_asv_df = pd.DataFrame(
            0,
            index=genome_ft.index,
            columns=["Estimate", "P-value", "BH-Corrected P-Value", "R2"]
        )

    # Loop through each ASV and merge with target_df
    for asv_name in genome_ft.index:
        #get each df in the results dictionary
        df = results[asv_name]
        
        #combine the empty_asv_df to that
        merged_pscores_df = empty_asv_df.combine_first(df)
        
        # Merge with target_df to include all sample columns
        merged_final_df = pd.concat([merged_pscores_df, target_df], axis=1)

        # Add the merged dataframe to the dictionary
        merged_list[asv_name] = merged_final_df

    return merged_list


def melt_correlation_results(results):
    """
    Melt a dictionary of correlation dataframes into a single dataframe.

    Parameters:
    results (dict): Dictionary where keys are ASV names, and values are dataframes 
                    with correlation results.

    Returns:
    pd.DataFrame: A single dataframe with all correlation results.
    """
    melted_results = []

    for asv_name, df in results.items():
        # Add Feature and Variable columns
        df = df.copy()
        df["Feature"] = df.index  # Index as Feature
        df["Variable"] = asv_name  # ASV name as Variable

        # Rearrange columns so Feature and Variable are first
        df = df[["Feature", "Variable"] + list(df.columns[:-2])]  # Reorder columns

        # Append to list
        melted_results.append(df)

    # Concatenate all dataframes
    final_df = pd.concat(melted_results, axis=0, ignore_index=True)

    return final_df

