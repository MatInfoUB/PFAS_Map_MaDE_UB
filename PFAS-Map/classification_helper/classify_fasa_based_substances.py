from pfas_object import PFAS
from .classify_pfaa_precursors import classify_fasas


def determine_alkyl_chain(smiles):
    c_set = set(list(smiles))
    return len(c_set) == 1 and 'C' in c_set


def classify_alkyl_fasas(pfas):  # N-alkyl perfluoroalkane sulfonamides
    smiles = pfas.get_rdkit_smiles()
    if 'NS(=O)(=O)' in smiles:
        split_smiles = smiles.split("NS(=O)(=O)", 1)   # Separate alkyl chain from PASF chains
        alkyl_chain = split_smiles[0]
        fasa_chain = "NS(=O)(=O)" + split_smiles[1]
        if determine_alkyl_chain(alkyl_chain) and classify_fasas(PFAS(fasa_chain)):
            return True
    return False


def classify_alkyl_fases(pfas):  # (N-alkyl) perfluoroalkane sulfonamidoethanols
    smiles = pfas.get_rdkit_smiles()
    if "N(CCO)S(=O)(=O)" in smiles:
        split_smiles = smiles.split("N(CCO)S(=O)(=O)", 1)
        alkyl_chain = split_smiles[0]
        fasa_chain = "NS(=O)(=O)" + split_smiles[1]
        if determine_alkyl_chain(alkyl_chain) and classify_fasas(PFAS(fasa_chain)):
            return True
    elif smiles.startswith("O=S(=O)(NCCO)"):
        split_smiles = smiles.split("O=S(=O)(NCCO)", 1)
        fasa_chain = "NS(=O)(=O)" + split_smiles[1]
        if classify_fasas(PFAS(fasa_chain)):
            return True
    return False


def classify_alkyl_fasacs(pfas):  # N-Alkyl perfluoroalkane sulfonamidoethyl acrylates
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("C=CC(=O)OCCN") and "S(=O)(=O)" in smiles:
        split_smiles = smiles.split("C=CC(=O)OCCN", 1)
        new_split_smiles = split_smiles[1].split("S(=O)(=O)", 1)
        alkyl_chain = new_split_smiles[0].strip("()")
        if alkyl_chain != new_split_smiles[0]:   # to prevent C=CC(=O)OCCNCCCS(=O)(=O)
            fasa_chain = "NS(=O)(=O)" + new_split_smiles[1]
            if determine_alkyl_chain(alkyl_chain) and classify_fasas(PFAS(fasa_chain)):
                return True
    return False


def classify_alkyl_fasmacs(pfas):  # N-Alkyl perfluoroalkane sulfonamidoethyl methacrylates
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("C=C(C)C(=O)OCCN") and "S(=O)(=O)" in smiles:
        split_smiles = smiles.split("C=C(C)C(=O)OCCN", 1)
        new_split_smiles = split_smiles[1].split("S(=O)(=O)", 1)
        alkyl_chain = new_split_smiles[0].strip("()")
        if alkyl_chain != new_split_smiles[0]:
            fasa_chain = "NS(=O)(=O)" + new_split_smiles[1]
            if determine_alkyl_chain(alkyl_chain) and classify_fasas(PFAS(fasa_chain)):
                return True
    return False


def classify_alkyl_fasaas(pfas):   # N-alkyl perfluoroalkane sulfonamido- acetic acids
    smiles = pfas.get_rdkit_smiles()
    if "N(CC(=O)O)S(=O)(=O)" in smiles:
        split_smiles = smiles.split("N(CC(=O)O)S(=O)(=O)", 1)
        alkyl_chain = split_smiles[0]
        fasa_chain = "NS(=O)(=O)" + split_smiles[1]
        if determine_alkyl_chain(alkyl_chain) and classify_fasas(PFAS(fasa_chain)):
            return True
    elif smiles.startswith("O=C(O)CNS(=O)(=O)"):
        split_smiles = smiles.split("O=C(O)CNS(=O)(=O)", 1)
        fasa_chain = "NS(=O)(=O)" + split_smiles[1]
        if classify_fasas(PFAS(fasa_chain)):
            return True
    return False


def classifying_fasa_precursors(pfas):
    if classify_alkyl_fasas(pfas):
        return ['FASA based PFAA precursors', 'N-Alkyl FASAs']
    if classify_alkyl_fases(pfas):
        return ['FASA based PFAA precursors', '(N-Alkyl) FASEs']
    if classify_alkyl_fasacs(pfas):
        return ['FASA based PFAA precursors', 'N-Alkyl FASACs']
    if classify_alkyl_fasmacs(pfas):
        return ['FASA based PFAA precursors', 'N-Alkyl FASMACs']
    if classify_alkyl_fasaas(pfas):
        return ['FASA based PFAA precursors', '(N-Alkyl) FASAAs']
    return None

