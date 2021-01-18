from pfas_object import PFAS
import pandas as pd
from .classify_pfaas import classifying_pfaas
from .classify_pfaa_precursors import classifying_pfaa_precursors
from .classify_other_perfluoroalkyls import classifying_other_perfluoroalkyls
from .classify_fasa_based_substances import classifying_fasa_precursors
from .classify_fluorotelomers import classifying_fluorotelomers
from rdkit import Chem


def replace_halogen(rdkit_smiles):
    smiles = rdkit_smiles.replace("Cl", "F")
    smiles = smiles.replace("Br", "F")
    smiles = smiles.replace("I", "F")
    return smiles


def generate_all_cf2_structure():
    atom_list = ['C', 'O', 'N', 'P', 'S']
    cf2_list = []
    for i in range(0, len(atom_list)):
        for j in range(i, len(atom_list)):
            cf2_list.append(atom_list[i] + 'C(F)(F)' + atom_list[j])
    return cf2_list


def generate_all_doublebond_structure():
    atom_list = ['C', 'O', 'N', 'P', 'S', 'F']
    doublebond_list = []
    for i in range(0, len(atom_list)):
        doublebond_list.append('C=C(F)' + atom_list[i])
        doublebond_list.append('O=C(F)' + atom_list[i])
    return doublebond_list


def determine_perfluoro_unit(rdkit_smiles):
    patt_1 = Chem.MolFromSmiles("C(F)(F)F")
    m = Chem.MolFromSmiles(rdkit_smiles)
    if m.HasSubstructMatch(patt_1):
        return True
    cf2_list = generate_all_cf2_structure()
    for cf2 in cf2_list:
        patt = Chem.MolFromSmiles(cf2)
        if m.HasSubstructMatch(patt):
            return True
    return False


def determine_doublebond_unit(rdkit_smiles):
    double_list = generate_all_doublebond_structure()
    m = Chem.MolFromSmiles(rdkit_smiles)
    for doublebond in double_list:
        patt = Chem.MolFromSmiles(doublebond)
        if m.HasSubstructMatch(patt):
            return True
    return False


def find_aromatic_fluorinated_carbon(rdkit_smiles):
    condition = "c(F)" in rdkit_smiles
    return condition


def classify_pfas_derivatives(rdkit_smiles):
    replaced_smiles = replace_halogen(rdkit_smiles)
    if determine_perfluoro_unit(replaced_smiles):
        return ['PFAS derivatives', 'PFAS halogen derivatives']
    if find_aromatic_fluorinated_carbon(rdkit_smiles):
        return ['PFAS derivatives', 'With fluorinated aromatic carbon']
    if determine_doublebond_unit(rdkit_smiles):
        return ['PFAS derivatives', 'With fluorinated C=C or C=O carbon']
    return ['Not PFAS', 'Not PFAS by current definition']


def find_pacf_pasf_substances(rdkit_smiles):
    pacf = Chem.MolFromSmiles("O=C(O)C(F)(F)C")
    pacf_2 = Chem.MolFromSmiles("O=C(O)C(F)(F)F")
    pasf = Chem.MolFromSmiles("O=S(=O)C(F)(F)C")
    pasf_2 = Chem.MolFromSmiles("O=S(=O)C(F)(F)F")

    m = Chem.MolFromSmiles(rdkit_smiles)
    if m.HasSubstructMatch(pacf) or m.HasSubstructMatch(pacf_2):
        return "PACF-based substances"
    if m.HasSubstructMatch(pasf) or m.HasSubstructMatch(pasf_2):
        return "PASF-based substances"
    return "Others"


def contain_silicon(rdkit_smiles):
    return "Si" in rdkit_smiles


def classify_silicons(rdkit_smiles):
    return ['Silicon PFASs', find_pacf_pasf_substances(rdkit_smiles)]


def determine_aromatics(rdkit_smiles):
    return "c" in rdkit_smiles


def classify_aromatics(rdkit_smiles):
    return ['Side-chain aromatics', find_pacf_pasf_substances(rdkit_smiles)]


def determine_cyclic(rdkit_smiles):
    for i in range(1, 6):
        if str(i) in rdkit_smiles:
            return True
    return False


def replace_cyclic(rdkit_smiles):
    for i in range(1, 6):
        rdkit_smiles = rdkit_smiles.replace(str(i), "(F)")
    return rdkit_smiles


def cyclic_class(classification):
    classification[0] = classification[0] + ", cyclic"
    classification[1] = classification[1] + ", cyclic"
    return classification


def classify_aliphatic(pfas):
    pfaas = classifying_pfaas(pfas)  # Classifying PFAAs
    if pfaas is not None:
        return pfaas

    precursors = classifying_pfaa_precursors(pfas)  # Classifying perfluoro PFAA precursors
    if precursors is not None:
        return precursors

    other_perfluoro = classifying_other_perfluoroalkyls(pfas)  # Classifying other perfluoroalkyls
    if other_perfluoro is not None:
        return other_perfluoro

    fasa = classifying_fasa_precursors(pfas)  # Classifying FASA based substances
    if fasa is not None:
        return fasa

    telomer = classifying_fluorotelomers(pfas)  # Classifying fluorotelomer substances
    if telomer is not None:
        return telomer
    return None


def classify_pfas_molecule(rdkit_smiles, exist=False, df_descriptors=pd.DataFrame()):
    print("Classifying", rdkit_smiles)

    if not determine_perfluoro_unit(rdkit_smiles):  # Classifying PFAS derivatives/non-PFASs
        return classify_pfas_derivatives(rdkit_smiles)

    if contain_silicon(rdkit_smiles):  # Classifying silicon PFASs
        return classify_silicons(rdkit_smiles)

    if determine_aromatics(rdkit_smiles):  # Classifying aromatics
        return classify_aromatics(rdkit_smiles)

    # Transform cyclic pfas to linear pfas and initiate PFAS class
    cylic = False
    if determine_cyclic(rdkit_smiles):
        cylic = True
        rdkit_smiles = replace_cyclic(rdkit_smiles)
        pfas = PFAS(rdkit_smiles)
        if len(pfas.get_descriptors()) == 0:
            return ['Unable to open ring(s)', 'Unable to open ring(s)']
    else:
        pfas = PFAS(rdkit_smiles, exist=exist, df_descriptors=df_descriptors)
        if len(pfas.get_descriptors()) == 0:
            return ['Complex structure', 'Complex structure']

    # Classification through traditional aliphatic PFAS classification system
    aliphatic_class = classify_aliphatic(pfas)
    if aliphatic_class is not None:
        if cylic:
            return cyclic_class(aliphatic_class)
        return aliphatic_class

    # For the rest of aliphatics, identify PACF/PASF derivatives
    derivatives_class = ['Other aliphatics', find_pacf_pasf_substances(rdkit_smiles)]
    if cylic:
        return cyclic_class(derivatives_class)
    return derivatives_class


if __name__ == "__main__":
    cf2_list = generate_all_cf2_structure()
    for cf2 in cf2_list:
        print(cf2)