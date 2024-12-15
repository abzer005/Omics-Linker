import streamlit as st
from src.common import *
import pandas as pd

from streamlit.components.v1 import html

page_setup()

st.title('Corromics')

c1, c2, c3 = st.columns([2, 2, 2])  # Adjust column ratios for centering
with c2:
    try:
        st.image("assets/corromics_icon.png", caption="Corromics Logo", use_container_width=True)
    except TypeError:
        st.image("assets/corromics_icon.png", caption="Corromics Logo", use_column_width=True)

# Introduction
st.markdown("""
<div style="color: #05577C; font-size: 25px;">
<b>Welcome to Corromics!</b> Before proceeding with your analysis, please take a moment to read this homepage.
</div>
""", unsafe_allow_html=True)

st.write(' ')
st.subheader('What is Corromics Used For?')

st.write("""
         To be filled""")

# Input Files
st.subheader('Input File Requirements')
st.write(""" 

         """)

st.write("""
### Data Preparation Essentials
To ensure smooth processing, follow these guidelines:

- Input files must include the `.mzML` extension in the feature quantification table and metadata table.
- **Metadata Table**:
  - Must include a column **filename**.
  - Can include attribute columns such as **replicates**, **sample type** (e.g., control, treatment), etc.
  - **Time**: A specific column for time is mandatory.
""")

st.markdown("""          
Example feature table:  
 
|feature_ID|sample1.mzML|sample2.mzML|blank.mzML|
|---|---|---|---|
|1|1000|1100|100|
|2|2000|2200|200|
""")

st.write(' ')

st.markdown("""        
Example meta data table:
            
|filename|Sample_Type|Time_Point|
|---|---|---|
|sample1.mzML|Sample|1h|
|sample2.mzML|Sample|2h|
|blank.mzML|Blank|N/A| 
""")

# Output Files
st.subheader('Output File Information')
st.write("""
     
Upon processing your data, ChemProp2 generates an output in the form of a CSV file. 
This file is an enhanced version of the node-pair information, now including ChemProp2 scores for each pair. 
The scores range from -1 to +1, providing a comprehensive score for every edge (or node pair) within the molecular network.
         
Key aspects of the output CSV include:
- **Score Range**: Each node pair gets a score between -1 and +1.
- **Score Interpretation**: 
  - The magnitude of the score indicates the strength of the potential transformation. 
  - The sign of the score (positive or negative) reveals the directionality of the transformation. 
  - For example, in a node pair A & B, the sign of the score will indicate whether the transformation is from A to B or vice versa.
""")

st.write("""
### Integration with Cytoscape
- You can download the results as a **zip file** containing the required **GraphML** and **style files** for Cytoscape visualization.
- To use the output in Cytoscape:
  1. Open the **GraphML** file in Cytoscape.
  2. Import the **style file**:
     - Navigate to **File > Import > Styles from File**.
     - Select and upload the `.xml` style file.
  3. Apply the imported style:
     - Go to the **Styles** panel in Cytoscape.  
     - By default, the style will be set to **default**.  
     - Use the drop-down menu to select the uploaded style (it may appear as **default_0** or another variation).  
     - The visualization will update based on the selected style.
        
         """)

# Subheader and Interactive Features
st.subheader('About the App Elements')
st.markdown("""
🔍 **How to Know If the App Is Running?**  
If you're performing a calculation or generating a figure and don't see any updates on the main page, 
check the **top-right corner** of the page. If you see the message **'RUNNING'**, the app is active and processing your request.  
            
💡 **All plots are interactive!**  
- Use your mouse to select areas and zoom in on specific regions.  
- Double-click on the plot to reset the zoom and return to the default view.  
- Save plots using the **camera icon** in the top-right corner of the plot. You can specify the image format (e.g., PNG) in the settings panel before saving.
""")

# Citation and Resources
st.subheader('Citation and Further Resources')
st.write('If you use Corromics in your research, please cite:')
st.markdown("""
            * [FBMN-STATS](https://fbmn-statsguide.gnps2.org/) - A statistical pipeline for downstream processing of FBMN results.
            * Pakkir Shah, A.K., Walter, A., Ottosson, F. et al. Statistical analysis of feature-based molecular networking results from non-targeted metabolomics data. Nat Protoc (2024). https://doi.org/10.1038/s41596-024-01046-3
            """
            )
            

# Add more links as needed

# Feedback Section
st.subheader("We Value Your Feedback")
st.markdown("""
            We welcome your feedback and suggestions to improve Corromics. Please feel free to create an issue on our GitHub repository to share your thoughts or report any issues you encounter. 
            Your input is invaluable in making the tool better for everyone.

            [Create an Issue on GitHub](https://github.com/abzer005/Omics-Linker/issues/new)
""")

# Contribution and Follow Us
st.subheader("Contribute and Follow Us")
st.markdown("""
- Interested in contributing? Check out the [GitHub page](https://github.com/abzer005).
- For more about our work, visit our [lab's GitHub page](https://github.com/Functional-Metabolomics-Lab).
""")

# Optional: Footer
st.markdown("---")
st.text("Corromics © 2024")