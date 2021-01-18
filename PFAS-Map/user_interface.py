import pandas as pd
import streamlit as st
from utils import train_pca_pc74, process_epa_pfas, plot_unknown_on_chemical_map, pca_for_unknown_pfas_pc74,\
    merge_epa_and_unknown, fit_transform_tsne, get_tsne_results_with_property, show_property_in_3d_tsne


@st.cache
def load_data():
    selected_descriptors = pd.read_csv("processed_data/List_of_selected_descriptors.csv")
    for_pca_training = pd.read_csv("processed_data/Precalculated_EPA_PFAS_descriptors_data.csv")\
        .set_index("RDKIT_SMILES")
    preload_epa_tsne_classes = pd.read_csv("processed_data/TSNE_results_with_classification_PC_74.csv")\
        .set_index("RDKIT_SMILES")
    preload_epa_tsne_classes.drop_duplicates(inplace=True)
    preload_epa_pca_classes = pd.read_csv("processed_data/Precalculated_EPA_PFAS_pca_classification_pc74.csv")\
        .set_index("RDKIT_SMILES")
    preload_epa_pca_classes.drop_duplicates(inplace=True)
    return selected_descriptors, for_pca_training, preload_epa_tsne_classes, preload_epa_pca_classes


@st.cache
def get_scaler_pca(for_pca_training):
    return train_pca_pc74(for_pca_training)


@st.cache
def run_tsne(df_for_tsne):
    return fit_transform_tsne(df_for_tsne)


st.sidebar.title("PFAS-Map")

st.sidebar.info("PFAS-Map is a database framework for the rapid screening of structure-function relationships in PFASs."
                " You can explore the classification of PFASs from US EPA PFAS Master List, classify an unknown PFAS, "
                "or explore your experimental or machine-learning data of PFASs through PFAS-Map.")
option = st.sidebar.selectbox("What would you like to do?", ["Explore the classification of PFASs",
                                                             "Classify a potential PFAS",
                                                             "Explore your PFAS property data"])
df_selected_descriptors, df_for_pca_training, df_preload_epa_tsne_classes, df_preload_epa_pca_classes = load_data()
scaler, pca = get_scaler_pca(df_for_pca_training)
df_epa_pca_normalized = process_epa_pfas(df_preload_epa_pca_classes.copy())


if option == "Explore the classification of PFASs":
    st.title("Explore the classification of PFASs")
    explore_option = st.selectbox("Would you like to explore all PFASs or a specific class of PFAS?",
                                       ["All", "A specific class"])
    df_epa_pfas_normalized = process_epa_pfas(df_preload_epa_tsne_classes.copy())
    if explore_option == "All":
        fig = plot_unknown_on_chemical_map(df_epa_pfas_normalized, None)
        st.plotly_chart(fig)

    if explore_option == "A specific class":
        class_list = ["PFAAs", "Non-PFAA perfluoroalkyl substances", "Perfluoroalkyl PFAA precursors",
                      "FASA-based PFAA precursors", "Fluorotelomer-based PFAA precursors",
                      "Side-chain fluorinated aromatic PFASs", "Other aliphatic PFASs", "PFASs containing Silicon",
                      "PFAS derivatives"]
        class_option = st.selectbox("Which class of PFASs you would like to explore?",
                                      class_list)
        # df_epa_pfas_normalized = process_epa_pfas(df_epa_pfas_tsne.copy())
        df_epa_selected_class = df_epa_pfas_normalized[df_epa_pfas_normalized['Class'] == class_option].copy()
        fig = plot_unknown_on_chemical_map(df_epa_selected_class, df_unknown_pfas=None, subclass=True)
        st.plotly_chart(fig)

if option == "Classify a potential PFAS":
    input_smiles = st.sidebar.text_input("Please input the SMILES of the PFAS you want to classify",
                                         "C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(C(C(C(F)(F)F)(F)F)(F)F)(F)F")
    df_smiles = pd.DataFrame(data=[input_smiles], columns=['SMILES'])
    df_unknown_pfas = pca_for_unknown_pfas_pc74(scaler, pca, df_smiles, df_selected_descriptors.copy())
    unknown_class = df_unknown_pfas.iloc[0, -2]
    unknown_subclass = df_unknown_pfas.iloc[0, -1]

    st.title("Classify a Potential PFAS")
    st.write("Classification results:")
    st.write("Class:", unknown_class)
    st.write("Subclass:", unknown_subclass)
    st.write("Plotting the potential PFAS in PFAS-Map...(it may take several minutes)")

    df_merged = merge_epa_and_unknown(df_epa_pca_normalized, df_unknown_pfas)
    df_tsne_transformed = run_tsne(df_merged)
    df_unknown_tsne = df_tsne_transformed[df_tsne_transformed.index.isin(list(df_unknown_pfas.index))]
    fig = plot_unknown_on_chemical_map(df_tsne_transformed, df_unknown_tsne)
    df_epa_unknown_class = df_tsne_transformed[df_tsne_transformed['Class'] == unknown_class].copy()
    fig_2 = plot_unknown_on_chemical_map(df_epa_unknown_class, df_unknown_tsne, subclass=True)
#
    plot_option = st.selectbox("Select the level of classification you would like to see", ["Classes", "Subclasses"])
    if plot_option == "Classes":
        st.plotly_chart(fig)
    if plot_option == "Subclasses":
        st.plotly_chart(fig_2)

if option == "Explore your PFAS property data":
    input_file = st.sidebar.file_uploader("Choose a csv file with PFAS property data, use 'NAME' as column name for "
                                          "compound's name, 'SMILES' as column name for compound's SMILES, "
                                          "and 'PROPERTY' as column name for compound's property data", type="csv")
    if input_file is not None:
        st.title("Explore your PFAS property data")
        df_input = pd.read_csv(input_file, usecols=["NAME", "SMILES", "PROPERTY"]).dropna(subset=['SMILES'])
        st.write("Classifying your PFASs and plotting their property data in PFAS-Map...(it may take several minutes)")
        df_join = get_tsne_results_with_property(df_input, df_for_pca_training, df_epa_pca_normalized,
                                                 df_selected_descriptors)
        fig_3D = show_property_in_3d_tsne(df_join, property_name='Property',
                                          hover_data=['PROPERTY', 'NAME', 'SMILES', 'Class', 'Subclass'])
        st.plotly_chart(fig_3D)