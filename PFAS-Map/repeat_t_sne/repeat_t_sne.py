import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot


def plot_classes_tsne(dfs, names, opacity, plotname):
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
            name=names[i]
        )
        trace_list.append(trace)

    layout = go.Layout(
        scene=dict(
            xaxis=dict(
                title=dict(
                    text='TSNE-PCA-1',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                # range=[0, 18],
            ),
            yaxis=dict(
                title=dict(
                    text='TSNE-PCA-2',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                # range=[0, 12],
            ),
            zaxis=dict(
                title=dict(
                    text='TSNE-PCA-3',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                # range=[-16, -8],
            ), ),
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
                size=20,
                color="black"
            ),
            bordercolor="Grey",
            borderwidth=1
        ),
        template="plotly_white"
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plot(fig, filename=plotname)


def plot_classes_tsne_2d(dfs, names, opacity, plotname, x=1, y=2):
    trace_list = []
    for i in range(len(dfs)):
        trace = go.Scatter(
            x=dfs[i].iloc[:, x-1],
            y=dfs[i].iloc[:, y-1],
            text=dfs[i].index,
            mode='markers',
            marker=dict(
                size=5,
                opacity=opacity[i]
            ),
            name=names[i]
        )
        trace_list.append(trace)

    layout = go.Layout(
        legend=go.layout.Legend(
            traceorder="normal",
            font=dict(
                family="sans-serif",
                size=20,
                color="black"
            ),
            bordercolor="Grey",
            borderwidth=1
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)

    fig.update_layout(template="plotly_white", font=dict(family="Arial", size=9),
        width=692, height=450, margin=dict(l=60, r=20, t=40, b=40), showlegend=False)
    fig.update_xaxes(automargin=True, zeroline=True, zerolinewidth=1, zerolinecolor='Black', mirror=True, ticks='outside',
                     showline=True, linecolor='Black',
                     # range=[-5, 20],
                     title=dict(
                    text='TSNE-PCA-'+str(x),
                    font=dict(family="sans-serif", size=12, color="black"))
                     )
    fig.update_yaxes(automargin=True, zeroline=True, zerolinewidth=1, zerolinecolor='Black', mirror=True, ticks='outside',
                     showline=True, linecolor='Black',
                     # range=[-20, 15],
                     title=dict(
                    text='TSNE-PCA-'+str(y),
                    font=dict(family="sans-serif", size=12, color="black"))
                     )
    fig.write_image(plotname)


def plot_main_classes(mainclass_traces):
    names = ["PFAAs", "Non-PFAA perfluoroalkyl substances", "Perfluoroalkyl PFAA precursors", "FASA-based PFAA precursors",
            "Fluorotelomer-based PFAA precursors", "Side-chain fluorinated aromatic PFASs", "Other aliphatic PFASs",
            "PFASs containing Silicon", "PFAS derivatives"]
    plot_classes_tsne(mainclass_traces, names, np.ones(len(names)), 't_sne_outputs/TSNE_Main_classes.html')
    # Plot 2d projection based on which two axes you want to see
    # plot_classes_tsne_2d(mainclass_traces, names, np.ones(len(names)),
    #                      't_sne_outputs/TSNE_Main_classes_tsne_1&2.pdf', x=1, y=2)


def plot_subclasses(class_trace, class_name, subclass_names):
    subclass_traces = []
    subclass_names = subclass_names
    for name in subclass_names:
        subclass_traces.append(class_trace[class_trace['Second_Class'].isin([name, name+', cyclic'])])
    plot_classes_tsne(subclass_traces, subclass_names, np.ones(len(subclass_names)),
                      't_sne_outputs/TSNE_'+class_name+'.html')
    # Plot 2d projection based on which two axes you want to see
    # plot_classes_tsne_2d(subclass_traces, subclass_names, np.ones(len(subclass_names)),
    #                      't_sne_outputs/TSNE_'+class_name+'tsne_2&3.pdf', x=2, y=3)


if __name__ == "__main__":
    # Step 1: Use PCA to reduce the number of components
    df_for_pca = pd.read_csv("../processed_data/Precalculated_EPA_PFAS_descriptors_data.csv").set_index("RDKIT_SMILES")
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df_for_pca)
    num_components = 74  # where explained variance reaches 70%
    pca = PCA(n_components=num_components)
    pca_result = pca.fit_transform(df_scaled)
    print('Cumulative explained variation for principal components: {}'.format(np.sum(pca.explained_variance_ratio_)))

    # # Step 2: Use t-SNE to visualize the PCA-reduced dataset in 3-dimensional space
    perp = 50   # hyper parameter for optimization
    step = 1000  # hyper parameter for optimization
    tsne = TSNE(n_components=3, verbose=0, perplexity=perp, n_iter=step, random_state=42)
    tsne_pca_result = tsne.fit_transform(pca_result)
    df_tsne_pca = pd.DataFrame(tsne_pca_result, columns=['TSNE-PCA-1','TSNE-PCA-2','TSNE-PCA-3'], index=df_for_pca.index)
    df_classification = pd.read_csv("../processed_data/Precalculated_EPA_PFAS_pca_classification.csv")\
        .set_index("RDKIT_SMILES")
    df_tsne_pca_result = df_tsne_pca.join(df_classification)
    tsne_result_output_name = "t_sne_outputs/TSNE_results_with_classification_PC_"+str(num_components)+".csv"
    df_tsne_pca_result.to_csv(tsne_result_output_name)
    # df_tsne_pca_result = pd.read_csv(tsne_result_output_name).set_index("RDKIT_SMILES")

    # Step 3: Make 3-dimensional plots for the main classes and subclasses
    trace_1 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['PFAAs', 'PFAAs, cyclic'])]
    trace_2 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Non-PFAA perfluoroalkyls',
                                                           'Non-PFAA perfluoroalkyls, cyclic'])]
    trace_3 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Perfluoro PFAA precursors',
                                                           'Perfluoro PFAA precursors, cyclic'])]
    trace_4 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['FASA based PFAA precursors',
                                                           'FASA based PFAA precursors, cyclic'])]
    trace_5 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Fluorotelomer PFAA precursors',
                                                           'Fluorotelomer PFAA precursors, cyclic'])]
    trace_6 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Side-chain aromatics'])]

    # Optional traces for not well-classified PFASs
    trace_7 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Other aliphatics', 'Other aliphatics, cyclic',
                                                           'Unable to open ring(s)'])]
    trace_8 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['Silicon PFASs'])]
    trace_9 = df_tsne_pca_result[df_tsne_pca_result['First_Class'].isin(['PFAS derivatives'])]

    # """Plot main classes on t-SNE"""
    # include only the well-classified PFASs
    # plot_main_classes([trace_1, trace_2, trace_3, trace_4, trace_5, trace_6])
    # include all PFASs classes
    plot_main_classes([trace_1, trace_2, trace_3, trace_4, trace_5, trace_6, trace_7, trace_8, trace_9])
    #
    # """Plot the subclasses of each main class on t-SNE"""
    plot_subclasses(trace_1, 'PFAAs', ['PFCAs', 'PFSAs', 'PFSiAs', 'PFECAs', 'PFESAs', 'PFPAs', 'PFPiAs'])
    plot_subclasses(trace_2, 'Non-PFAA perfluoroalkyl substances', ['PFAs', 'PFAenes', 'PFACs', 'PFAKs'])
    plot_subclasses(trace_3, 'Perfluoroalkyl PFAA precursors', ['PASFs', 'FASAs', 'PAFs', 'PFAIs', 'PFALs'])
    plot_subclasses(trace_4, 'FASA-based PFAA precursors', ['N-Alkyl FASAs', '(N-Alkyl) FASEs', 'N-Alkyl FASACs',
                                                            'N-Alkyl FASMACs', '(N-Alkyl) FASAAs'])
    plot_subclasses(trace_5, "Fluorotelomer-based PFAA precursors", ['n:2 FTIs', 'n:2 FTOs', 'n:2 FTOHs', 'n:2 FTACs',
                                                                     'n:2 FTMACs', 'n:2 monoPAPs', 'n:2 diPAPs',
                                                                     'n:2 FTALs', 'n:2 FTUALs', 'n:2 FTCAs',
                                                                     'n:2 FTUCAs', 'n:2 FTSAs', 'SFAs', 'SFAenes',
                                                                     'n:3 Acids', 'n:3 UAcids', 'n:1 FTOHs'])
    # n:2 fluorotelomers only
    # plot_subclasses(trace_5, "Fluorotelomer-based PFAA precursors", ['n:2 FTIs', 'n:2 FTOs', 'n:2 FTOHs', 'n:2 FTACs',
    #                                                                  'n:2 FTMACs', 'n:2 monoPAPs', 'n:2 diPAPs',
    #                                                                  'n:2 FTALs', 'n:2 FTUALs', 'n:2 FTCAs',
    #                                                                  'n:2 FTUCAs', 'n:2 FTSAs'])
    plot_subclasses(trace_6, 'Side-chain fluorinated aromatic PFASs', ['PACF-based substances', 'PASF-based substances'
                    , 'Others'])
