# Import necessary libraries
import streamlit as st
import pandas as pd
from src.common import *        # Importing common functionalities
from src.fileselection import * # Importing file selection functionalities
from src.correlation import *
from src.fdr import *
import random
from multiprocessing import Pool
from datetime import datetime

st.markdown("### Metabolomics-Metagenomics Data Combined")

#st.session_state['metabolome_md']
#st.session_state['genomics_md']

if 'metabolome_ft' in st.session_state and 'genomics_ft' in st.session_state:
    met_ft = st.session_state['metabolome_ft']
    gen_ft = st.session_state['genomics_ft']
    gen_ft = gen_ft[~(gen_ft.eq(0).all(axis=1))]

    # Combine the DataFrames
    target_df = combine_dataframes(met_ft, gen_ft)

    #Create a Decoy set
    random.seed(42)

    # Permute the values of genomics_ft
    decoy_gen_df = gen_ft.apply(np.random.permutation)
    decoy_gen_df += np.random.normal(0, 0.1, decoy_gen_df.shape)  # Add noise
    decoy_df = combine_dataframes(met_ft, decoy_gen_df)

    # Shuffle rows and columns of genomics_ft
    #shuffled_columns = random.sample(list(gen_ft.columns), len(gen_ft.columns))
    #shuffled_rows = random.sample(list(gen_ft.index), len(gen_ft.index))

    #gen_ft_decoy = gen_ft.loc[shuffled_rows, shuffled_columns]
    #decoy_df = combine_dataframes(met_ft, gen_ft_decoy)

    # Display the combined DataFrame in Streamlit
    with st.expander(f"Target Dataframe {target_df.shape}"):
                st.dataframe(target_df)
    with st.expander(f"Decoy Dataframe {decoy_df.shape}"):
                st.dataframe(decoy_df)

    # Perform Correlation ########################################################
    # Record the start time
    # start_time = datetime.now()
    # st.write(f"Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    #Give the time message to the user
    estimate_run_time(met_ft, gen_ft)
    
    with st.spinner("Calculating correlations..."):
        target_results = calculate_correlations_parallel(target_df, met_ft, gen_ft)
        decoy_results = calculate_correlations_parallel(decoy_df, met_ft, gen_ft)
        #merged_list = merge_asv_correlation_results(results, target_df, gen_ft

    melted_target = melt_correlation_results(target_results)
    melted_decoy = melt_correlation_results(decoy_results)

    st.session_state['Target_scores'] = melted_target
    st.session_state['Decoy_scores'] = melted_decoy

    with st.expander(f"Correlation Scores of Target Dataframe {melted_target.shape}"):
                st.dataframe(melted_target)
    with st.expander(f"Correlation Scores of Decoy Dataframe {melted_decoy.shape}"):
                st.dataframe(melted_decoy)
    
    st.markdown('### False Discovery Rate')

    if 'Target_scores' in st.session_state and 'Decoy_scores' in st.session_state:
        
        target_scores = st.session_state['Target_scores']
        decoy_scores = st.session_state['Decoy_scores']

        st.write('Select the positive and negative cutoffs for the correlation scores based on your FDR-curve')
        overall_fdr_table, fig_histogram, fig_fdr = calculate_fdr(target_scores, decoy_scores, score_range=(-1, 1), bin_size=0.001)

        st.plotly_chart(fig_histogram)
        st.plotly_chart(fig_fdr)

        # # Record the start time
        # end_time = datetime.now()
        # st.write(f"Start Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")




else:
    st.warning("Please input the data in the first page to continue the analysis here")
