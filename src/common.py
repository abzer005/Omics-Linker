import streamlit as st
import pandas as pd
import io
import uuid

def clear_cache_button():
   if st.button("Clear Cache"):
        # Clear cache for both newer and older Streamlit versions
        if hasattr(st, "cache_data"):
            st.cache_data.clear()
        if hasattr(st, "cache_resource"):
            st.cache_resource.clear()
        st.success("Cache cleared!")

    # initialize global session state variables if not already present
    # DataFrames

def v_space(n, col=None):
    for _ in range(n):
        if col:
            col.write("")
        else:
            st.write("")

def page_setup():
    # streamlit configs
    st.set_page_config(
        page_title="Corromics",
        page_icon="assets/corromics_icon.png",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )
    for key in dataframe_names:
        if key not in st.session_state:
            st.session_state[key] = pd.DataFrame()
    if "data_preparation_done" not in st.session_state:
        st.session_state["data_preparation_done"] = False

    with st.sidebar:
        with st.expander("⚙️ Settings", expanded=True):
            st.selectbox(
                "image export format",
                ["svg", "png", "jpeg", "webp"],
                key="image_format",
            )
        v_space(1)
        # Add the clear cache button
        clear_cache_button()
        v_space(1)
        
        try:
            st.image("https://raw.githubusercontent.com/abzer005/Omics-Linker/main/streamlit/assets/corromics_icon.png",
                     use_container_width=True)
        except TypeError:
            st.image("https://raw.githubusercontent.com/abzer005/Omics-Linker/main/streamlit/assets/corromics_icon.png", 
                     use_column_width=True)

dataframe_names = ("md",
                   "ft",
                   "nw",
                   "an_gnps",
                   "an_analog")


def reset_dataframes():
    for key in dataframe_names:
        st.session_state[key] = pd.DataFrame()



def open_df(file):
    separators = {"txt": "\t", "tsv": "\t", "csv": ","}
    try:
        if type(file) == str:
            ext = file.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        else:
            ext = file.name.split(".")[-1]
            if ext != "xlsx":
                df = pd.read_csv(file, sep=separators[ext])
            else:
                df = pd.read_excel(file)
        
        # sometimes dataframes get saved with unnamed index, that needs to be removed
        if "Unnamed: 0" in df.columns:
            df.drop("Unnamed: 0", inplace=True, axis=1)
        return df
    except:
        return pd.DataFrame()

def show_table(df, title="", col="", download=True):
    if col:
        col = col
    else:
        col = st
    if download:
        col.download_button(
            f"Download Table",
            df.to_csv(sep="\t").encode("utf-8"),
            title.replace(" ", "-") + ".tsv",
            key=uuid.uuid1(),
        )
    col.dataframe(df, use_container_width=True)


def show_fig(fig, download_name, container_width=True):

    # Set default image format to 'svg' if not specified in session state
    image_format = st.session_state.get('image_format', 'svg')

    st.plotly_chart(
        fig,
        use_container_width=container_width,
        config={
            "displaylogo": False,
            "modeBarButtonsToRemove": [
                "zoom",
                "pan",
                "select",
                "lasso",
                "zoomin",
                "autoscale",
                "zoomout",
                "resetscale",
            ],
            "toImageButtonOptions": {
                "filename": download_name,
                "format": image_format,
            },
        },
    )


def download_plotly_figure(fig, filename="", col=""):
    buffer = io.BytesIO()
    fig.write_image(file=buffer, format="png")

    if col:
        col.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )
    else:
        st.download_button(
            label=f"Download Figure",
            data=buffer,
            file_name=filename,
            mime="application/png",
        )
