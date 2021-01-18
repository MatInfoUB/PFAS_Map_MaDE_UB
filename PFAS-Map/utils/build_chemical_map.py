import pandas as pd
from rdkit_helper import rdkit_smiles_from_input_smiles
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from padel_helper import calculate_padel_descriptors_fingerprints
import plotly.express as px
import plotly.graph_objs as go
from classification_helper import classify_pfas_molecule


def normalize_classes(initial_class):
    if initial_class in ['PFAAs', 'PFAAs, cyclic']:
        return "PFAAs"
    if initial_class in ['Non-PFAA perfluoroalkyls', 'Non-PFAA perfluoroalkyls, cyclic']:
        return "Non-PFAA perfluoroalkyl substances"
    if initial_class in ['Perfluoro PFAA precursors', 'Perfluoro PFAA precursors, cyclic']:
        return "Perfluoroalkyl PFAA precursors"
    if initial_class in ['FASA based PFAA precursors', 'FASA based PFAA precursors, cyclic']:
        return "FASA-based PFAA precursors"
    if initial_class in ['Fluorotelomer PFAA precursors', 'Fluorotelomer PFAA precursors, cyclic']:
        return "Fluorotelomer-based PFAA precursors"
    if initial_class in ['Side-chain aromatics']:
        return "Side-chain fluorinated aromatic PFASs"
    if initial_class in ['Other aliphatics', 'Other aliphatics, cyclic', 'Unable to open ring(s)']:
        return "Other aliphatic PFASs"
    if initial_class in ['Silicon PFASs']:
        return "PFASs containing Silicon"
    if initial_class in ['PFAS derivatives']:
        return "PFAS derivatives"
    if initial_class in ['Not PFAS']:
        return "Not PFAS"
    return initial_class


def normalize_subclass(subclass):
    return subclass.replace(", cyclic", "")


def process_epa_pfas(df_epa_pfas):
    df_epa_pfas['Class'] = df_epa_pfas['First_Class'].map(lambda x: normalize_classes(x))
    df_epa_pfas['Subclass'] = df_epa_pfas['Second_Class'].map(lambda x: normalize_subclass(x))
    df_epa_pfas.drop(columns=['First_Class', 'Second_Class'], inplace=True)
    return df_epa_pfas


def get_rdkit_smiles(df_input):
    df_input['RDKIT_SMILES'] = df_input['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
    df = df_input.dropna(subset=['RDKIT_SMILES'])
    return df


def get_descriptors_data(df_smiles):
    rdkit_smiles_list = list(df_smiles['RDKIT_SMILES'])
    list_len = len(rdkit_smiles_list)
    df_descriptors_all = pd.DataFrame()
    for i in range(0, list_len):
        rdkit_smiles = rdkit_smiles_list[i]
        descriptors = calculate_padel_descriptors_fingerprints(rdkit_smiles)
        if len(descriptors) != 0:
            df_descriptors = pd.DataFrame.from_dict(descriptors, orient='index').T
            df_descriptors.insert(0, column='RDKIT_SMILES', value=rdkit_smiles)
            df_descriptors_all = df_descriptors_all.append(df_descriptors)
    return df_descriptors_all


def train_pca_pc74(df_for_pca_train):
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df_for_pca_train)
    pca = PCA(n_components=74)
    pca.fit(df_scaled)
    return scaler, pca


def transform_pca_pc74(scaler, pca, df_selected_descriptors, df_descriptors_all):
    selected_descriptors = list(df_selected_descriptors['Descriptors'])
    df_descriptors = df_descriptors_all[selected_descriptors].copy()
    df_scaled = scaler.transform(df_descriptors)
    pca_transformed = pca.transform(df_scaled)
    columns = []
    for i in range(74):
        columns.append('PC-'+str(i+1))
    df_pca_transformed = pd.DataFrame(pca_transformed, columns=columns,
                                      index=df_descriptors_all['RDKIT_SMILES'])
    return df_pca_transformed


def pca_for_unknown_pfas_pc74(scaler, pca, df_smiles, df_selected_descriptors):
    df_rdkit_smiles = get_rdkit_smiles(df_smiles)
    df_descriptors_all = get_descriptors_data(df_rdkit_smiles)
    transformed_result = transform_pca_pc74(scaler, pca, df_selected_descriptors, df_descriptors_all)
    df_unknown_pfas = transformed_result.reset_index().copy()
    df_unknown_pfas['Classification'] = df_unknown_pfas['RDKIT_SMILES'].map(
        lambda x: classify_pfas_molecule(x, True, df_descriptors_all))
    df_unknown_pfas['Class'] = df_unknown_pfas['Classification'].map(lambda x: normalize_classes(x[0]))
    df_unknown_pfas['Subclass'] = df_unknown_pfas['Classification'].map(lambda x: normalize_subclass(x[1]))
    df_unknown_pfas.drop(columns=['Classification'], inplace=True)
    df_unknown_pfas.set_index("RDKIT_SMILES", inplace=True)
    return df_unknown_pfas


def merge_epa_and_unknown(df_epa_precalculated, df_unknown_pfas):
    df_unknown_pfas_distinct = df_unknown_pfas[~df_unknown_pfas.index.isin(list(df_epa_precalculated.index))]
    df_merged = pd.concat([df_epa_precalculated, df_unknown_pfas_distinct])
    return df_merged


def fit_transform_tsne(df_merged):
    tsne = TSNE(n_components=3, verbose=0, perplexity=50, n_iter=1000, random_state=42)
    df_pca_train_and_test = df_merged.iloc[:, 0:74]
    tsne_transformed = tsne.fit_transform(df_pca_train_and_test)
    df_tsne_transformed = pd.DataFrame(tsne_transformed, columns=['TSNE-PCA-1', 'TSNE-PCA-2', 'TSNE-PCA-3'],
                                      index=df_pca_train_and_test.index)
    df_tsne_transformed['Class'] = df_merged['Class']
    df_tsne_transformed['Subclass'] = df_merged['Subclass']
    return df_tsne_transformed


def plot_unknown_on_chemical_map(df_epa_pfas, df_unknown_pfas, subclass=False):
    dfs = []
    opacity = []
    if not subclass:
        names = df_epa_pfas['Class'].unique().tolist()
        if "Not PFAS" in names:
            names.remove("Not PFAS")
        for name in names:
            dfs.append(df_epa_pfas[df_epa_pfas['Class'] == name].copy())
            if name in ['Other aliphatic PFASs', 'PFASs containing Silicon', 'PFAS derivatives']:
                opacity.append(0.1)
            else:
                opacity.append(1)
    else:
        names = df_epa_pfas['Subclass'].unique().tolist()
        for name in names:
            dfs.append(df_epa_pfas[df_epa_pfas['Subclass'] == name].copy())
            opacity.append(1)
    trace_list = []

    for i in range(len(dfs)):
        trace = go.Scatter3d(
            x=dfs[i].iloc[:, 0],
            y=dfs[i].iloc[:, 1],
            z=dfs[i].iloc[:, 2],
            text=dfs[i].index,
            mode='markers',
            marker=dict(
                size=5,
                opacity=opacity[i]
            ),
            name=names[i],
            hoverlabel=dict(namelength=-1)
        )
        trace_list.append(trace)
    if df_unknown_pfas is not None:
        trace_unknown = go.Scatter3d(
            x=df_unknown_pfas.iloc[:, 0],
            y=df_unknown_pfas.iloc[:, 1],
            z=df_unknown_pfas.iloc[:, 2],
            text=df_unknown_pfas.index,
            mode='markers',
            marker=dict(
                size=8,
                opacity=1,
                color='black'
            ),
            name="The PFAS you input",
            hoverlabel=dict(namelength=-1)
        )
        trace_list.append(trace_unknown)

    layout = go.Layout(
        scene=dict(
            xaxis=dict(
                title=dict(
                    # text='PC1',
                    text='TSNE-PCA-1',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                # range=[-25, 50]
            ),
            yaxis=dict(
                title=dict(
                    # text='PC2',
                    text='TSNE-PCA-2',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
            ),
            zaxis=dict(
                title=dict(
                    # text='PC3',
                    text='TSNE-PCA-3',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                # range=[-12, 10]
            )
        ),
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        legend=go.layout.Legend(
            traceorder="normal",
            font=dict(
                family="sans-serif",
                size=12,
                color="black"
            ),
            bordercolor="Grey",
            borderwidth=1
        ),
        template="plotly_white"
    )

    fig = go.Figure(data=trace_list, layout=layout)
    fig.update_layout(legend_orientation="h")
    return fig


def get_tsne_results_with_property(df_input, df_for_pca, df_epa_precalculated, df_selected_descriptors):
    df_rdkit_smiles = get_rdkit_smiles(df_input)
    scaler, pca = train_pca_pc74(df_for_pca)
    df_pca_transformed = pca_for_unknown_pfas_pc74(scaler, pca, df_input, df_selected_descriptors)
    df_merged = merge_epa_and_unknown(df_epa_precalculated, df_pca_transformed)
    df_tsne = fit_transform_tsne(df_merged)
    df_join = df_rdkit_smiles.set_index("RDKIT_SMILES").join(df_tsne)
    return df_join


def show_property_in_3d_tsne(df_join, property_name, hover_data):
    fig = px.scatter_3d(df_join, x="TSNE-PCA-1", y="TSNE-PCA-2", z="TSNE-PCA-3", hover_data=hover_data,
                        color='PROPERTY', color_continuous_scale='bluered')
    fig.update_traces( marker_size=6)
    fig.update_layout(coloraxis_colorbar=dict(
        title=property_name), template="plotly_white"
        # , scene=dict(zaxis=dict(range=[-22, 0]))
    )
    return fig


def show_property_in_2d_tsne(df_join, property_name, hover_data, x="TSNE-PCA-1", y="TSNE-PCA-2" ):
    fig = px.scatter(df_join, x=x, y=y, hover_data=hover_data,
                     # text='NAME',
                     color='PROPERTY', color_continuous_scale='bluered')
    fig.update_traces(marker_size=9)
    # fig.update_traces(mode='markers+text', textposition="top center", marker_size=12)
    fig.update_layout(template="plotly_white", font=dict(family="Arial", size=9),
        width=692, height=450, margin=dict(l=60, r=20, t=40, b=40), showlegend=False)
    fig.update_layout(coloraxis_colorbar=dict(
        title=property_name), template="plotly_white")
    fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Black', mirror=True, ticks='outside',
                     showline=True, linecolor='Black')
    fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Black', mirror=True, ticks='outside',
                     showline=True, linecolor='Black')
    return fig


# Old functions for PC=3
# def train_pca(df_for_pca):
#     scaler = StandardScaler()
#     df_scaled = scaler.fit_transform(df_for_pca)
#     pca = PCA(n_components=3)
#     pca.fit(df_scaled)
#     return scaler, pca
#
#
# def transform_pca(scaler, pca, df_selected_descriptors, df_descriptors_all):
#     selected_descriptors = list(df_selected_descriptors['Descriptors'])
#     df_descriptors = df_descriptors_all[selected_descriptors].copy()
#     df_scaled = scaler.transform(df_descriptors)
#     pca_transformed = pca.transform(df_scaled)
#     df_pca_transformed = pd.DataFrame(pca_transformed, columns=['PC-1', 'PC-2', 'PC-3'],
#                                       index=df_descriptors_all['RDKIT_SMILES'])
#     return df_pca_transformed
#
#
# def pca_for_unknown_pfas(scaler, pca, df_smiles, df_selected_descriptors):
#     df_rdkit_smiles = get_rdkit_smiles(df_smiles)
#     df_descriptors_all = get_descriptors_data(df_rdkit_smiles)
#     transformed_result = transform_pca(scaler, pca, df_selected_descriptors, df_descriptors_all)
#     df_unknown_pfas = transformed_result.reset_index().copy()
#     df_unknown_pfas['Classification'] = df_unknown_pfas['RDKIT_SMILES'].map(
#         lambda x: classify_pfas_molecule(x, True, df_descriptors_all))
#     df_unknown_pfas['Class'] = df_unknown_pfas['Classification'].map(lambda x: normalize_classes(x[0]))
#     df_unknown_pfas['Subclass'] = df_unknown_pfas['Classification'].map(lambda x: normalize_subclass(x[1]))
#     df_unknown_pfas.drop(columns=['Classification'], inplace=True)
#     df_unknown_pfas.set_index("RDKIT_SMILES", inplace=True)
#     return df_unknown_pfas
#
#
# def show_property_in_3d(df_join, property_name, hover_data):
#     fig = px.scatter_3d(df_join, x="PC-1", y="PC-2", z="PC-3", hover_data=hover_data,
#                         color='PROPERTY', color_continuous_scale='bluered')
#     fig.update_traces( marker_size=6)
#     fig.update_layout(coloraxis_colorbar=dict(
#         title=property_name), template="plotly_white")
#     return fig
#
#
# def show_property_in_2d(df_join, property_name, hover_data):
#     fig = px.scatter(df_join, x="PC-1", y="PC-3", hover_data=hover_data,
#                         color='PROPERTY', color_continuous_scale='bluered')
#     fig.update_traces(marker_size=12)
#     fig.update_layout(coloraxis_colorbar=dict(
#         title=property_name), template="plotly_white")
#     fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Black', mirror=True, ticks='outside',
#                      showline=True, linecolor='Black')
#     fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Black', mirror=True, ticks='outside',
#                      showline=True, linecolor='Black')
#     return fig
#
#
# def get_pca_results_with_property(df_input, df_for_pca, df_selected_descriptors):
#     df_rdkit_smiles = get_rdkit_smiles(df_input)
#     scaler, pca = train_pca(df_for_pca)
#     df_pca_transformed = pca_for_unknown_pfas(scaler, pca, df_input, df_selected_descriptors)
#     df_join = df_rdkit_smiles.set_index("RDKIT_SMILES").join(df_pca_transformed)
#     return df_join