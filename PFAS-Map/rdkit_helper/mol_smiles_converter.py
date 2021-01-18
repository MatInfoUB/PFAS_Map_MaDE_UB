from rdkit.Chem import SaltRemover, MolToSmiles, MolFromSmiles, RemoveHs
import numpy as np
from .Charge_neutralizer import NeutraliseCharges
from. isotope_remover import MolWithoutIsotopesToSmiles

SALTS_FILE = "salts/Salts_from_epa_list.txt"


def mol_from_smiles(smiles):
    try:
        mol = MolFromSmiles(smiles)
    except BaseException as error:
        print('An exception occurred: {}'.format(error))
        mol = None
    if mol is None:
        return None
    return MolWithoutIsotopesToSmiles(mol)


def rdkit_smiles_from_mol(mol):
    res = SaltRemover.SaltRemover(defnFilename=SALTS_FILE).StripMol(mol)
    smiles = MolToSmiles(res)
    smiles, _ = NeutraliseCharges(smiles)
    if "." in smiles:  # Check if structure after salt removal is chelate
        split_smiles = smiles.split(".", 5)
        if len(set(split_smiles)) == 1:
            return split_smiles[0]
    return smiles


def rdkit_smiles_from_input_smiles(smiles):
    mol = mol_from_smiles(smiles)
    if mol is not None:
        return rdkit_smiles_from_mol(mol)
    return np.nan

