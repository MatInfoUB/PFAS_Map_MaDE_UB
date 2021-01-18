import pandas as pd
import numpy as np
from rdkit_helper import rdkit_smiles_from_input_smiles
from padel_helper import calculate_padel_descriptors_fingerprints
from classification_helper import classify_pfas_molecule
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def prepare_rdkit_smiles(input_smiles_file):
    df = pd.read_excel(input_smiles_file, usecols=['SMILES'], na_values=['-', ''])
    # STEP ONE Remove PFASs with NO SMILES
    df.dropna(inplace=True)
    # STEP TWO Convert Input SMILES to RDKit SMILES (including removing salts and extract PFAS units from complexes)
    df['RDKIT_SMILES'] = df['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
    # STEP THREE Remove PFASs with INVALID SMILES (CANNOT CONVERT TO VALID RDKIT SMILES)
    df.dropna(axis=0, inplace=True)
    # STEP FOUR Remove Duplicates
    df.drop_duplicates(subset=['RDKIT_SMILES'], inplace=True)
    # STEP FIVE Remove Complexes containing more than one type of PFAS units
    df = df[~df['RDKIT_SMILES'].str.contains("\.")]
    # STEP Six Remove Isomeric SMILES (To simplify classification)
    df = df[~df['SMILES'].str.contains('@')]  # Remove isomeric smiles
    return df


if __name__ == "__main__":
    # create rdkit smiles
    smiles_csv = "input_data/EPA_PFAS_MASTERLIST-2020-05-04.xlsx"
    df_smiles = prepare_rdkit_smiles(smiles_csv)

    # calculate padel descriptors
    smiles_list = list(df_smiles['SMILES'])
    rdkit_smiles_list = list(df_smiles['RDKIT_SMILES'])
    list_len = len(smiles_list)
    df_descriptors_all = pd.DataFrame()
    df_timeouts_all = pd.DataFrame()
    for i in range(0, list_len):
        smiles = smiles_list[i]
        rdkit_smiles = rdkit_smiles_list[i]
        print("Calculating Descriptors for", rdkit_smiles)
        descriptors = calculate_padel_descriptors_fingerprints(rdkit_smiles)
        if len(descriptors) != 0:
            df_descriptors = pd.DataFrame.from_dict(descriptors, orient='index').T
            df_descriptors.insert(0, column='RDKIT_SMILES', value=rdkit_smiles)
            df_descriptors.insert(0, column='SMILES', value=smiles)
            df_descriptors_all = df_descriptors_all.append(df_descriptors)
        else:
            df_timeouts = pd.DataFrame()
            df_timeouts = df_timeouts.append({'SMILES': smiles, 'RDKIT_SMILES': rdkit_smiles}, ignore_index=True)
            df_timeouts_all = df_timeouts_all.append(df_timeouts)
            print("Warning: the SMILES causes PADELPY running time error")
        list_len = list_len - 1
        print(str(list_len) + " PFASs left")

    # process descriptors data
    df = df_descriptors_all.replace(["Infinity", "-Infinity", ""], np.nan)
    column_list = list(df.columns)
    index_3d_begin = column_list.index("TDB1u")
    index_3d_end = column_list.index("Ds")
    column_3d_list = column_list[index_3d_begin:index_3d_end + 1].copy()
    df = df.drop(columns=column_3d_list)
    df_column_processed = df.dropna(axis=1, thresh=df.shape[0] - 2)
    df_row_processed = df_column_processed.dropna(axis=0)

    df_for_pca = df_row_processed.drop(columns=['SMILES']).set_index('RDKIT_SMILES')
    df_screened_descriptors = pd.DataFrame(list(df_for_pca.columns), columns=['Descriptors'])
    df_screened_descriptors.to_csv("output_data/List_of_selected_descriptors_repeat.csv", index=False)
    df_for_pca.to_csv("output_data/Precalculated_EPA_PFAS_descriptors_data_repeat.csv")

    # classify all epa pfass
    df = df_for_pca.copy()
    df.reset_index(inplace=True)
    df_copy = df.copy()
    df['Classification'] = df['RDKIT_SMILES'].map(lambda x: classify_pfas_molecule(x, True, df_copy))
    df_classification = df[['RDKIT_SMILES', 'Classification']].copy()
    df_classification['First_Class'] = df_classification['Classification'].map(lambda x: x[0])
    df_classification['Second_Class'] = df_classification['Classification'].map(lambda x: x[1])
    df_classification.drop(columns='Classification', inplace=True)
    df_classification.to_csv("output_data/EPA_PFAS_CLASSIFICATION_repeat.csv", index=False)

    # temporary code for pca 74
    # df_for_pca = pd.read_csv("output_data/Precalculated_EPA_PFAS_descriptors_data.csv").set_index("RDKIT_SMILES")
    # df_classification = pd.read_csv("output_data/0506_CLASSIFICATION_2020_Refactored.csv")

    # get pca result
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df_for_pca)
    # pca = PCA(n_components=3)
    pca = PCA(n_components=74)
    pca_result = pca.fit_transform(df_scaled)
    columns = []
    for i in range(74):
        columns.append('PC-'+str(i+1))
    # df_pca_result = pd.DataFrame(pca_result, columns=['PC-1', 'PC-2', 'PC-3'], index=df_for_pca.index)
    df_pca_result = pd.DataFrame(pca_result, columns=columns, index=df_for_pca.index)
    df_join = df_pca_result.join(df_classification.set_index("RDKIT_SMILES")).reset_index()
    df_join.to_csv("output_data/Precalculated_EPA_PFAS_pca_classification_pc74_repeat.csv", index=False)


