from padelpy import from_smiles
from padelpy import padeldescriptor


def generate_padel_descriptors(smiles):
    return from_smiles(smiles)


def fetch_padel_descriptors(smiles, df):
    df_descriptors = df[df['RDKIT_SMILES'] == smiles].copy()
    # df_descriptors.drop(columns=['SMILES', 'RDKIT_SMILES'], inplace=True)
    df_descriptors.drop(columns=['RDKIT_SMILES'], inplace=True)
    return df_descriptors.to_dict('records')[0]


def calculate_padel_descriptors_fingerprints(smiles):
    # padeldescriptor(removesalt=True, detectaromaticity=True, d_3d=False, fingerprints=True, standardizetautomers=True,
    #                 standardizenitro=True, maxruntime=30000)
    try:
        descriptors = from_smiles(smiles, descriptors=True, fingerprints=True)
    except RuntimeError:
        descriptors = {}
    return descriptors

