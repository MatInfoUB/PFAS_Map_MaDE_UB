from .classify_pfaa_precursors import classify_pfais
from .classify_fasa_based_substances import determine_alkyl_chain
from pfas_object import PFAS


# Fluorotelomer PFAA precursors
def classify_ftis(pfas):  # n:2 Fluorotelomer iodides
    smiles = pfas.get_rdkit_smiles()
    if "CCI" in smiles:
        smiles = smiles.replace("CCI", "I")   # Remove the CH2CH2
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftos(pfas):  # n:2 Fluorotelomer olefins
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("C=CC"):
        smiles = smiles.replace("C=C", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftohs(pfas):  # n:2 Fluorotelomer alcohols
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("OCCC"):
        smiles = smiles.replace("OCC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftacs(pfas):  # n:2 Fluorotelomer acrylates
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("C=CC(=O)OCCC"):
        smiles = smiles.replace("C=CC(=O)OCC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftmacs(pfas):  # n:2 Fluorotelomer methacrylates
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("C=C(C)C(=O)OCCC"):
        smiles = smiles.replace("C=C(C)C(=O)OCC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_monoesters(pfas):  # n:2 Polyfluoroalkyl phosphoric acid esters, monoester
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=P(O)(O)OCCC"):
        smiles = smiles.replace("O=P(O)(O)OCC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_diesters(pfas):  # n:2 Polyfluoroalkyl phosphoric acid esters, diester
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=P(O)(OCCC"):
        smiles = smiles.replace("O=P(O)(OCC", "")
        split_smiles = smiles.split(")OCC", 1)
        if len(set(split_smiles)) == 1:
            new_smiles = "I" + split_smiles[0]
            if classify_pfais(PFAS(new_smiles)):
                return True
    return False


def classify_ftals(pfas):  # n:2 Fluorotelomer aldehydes
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=CCC"):
        smiles = smiles.replace("O=CC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftuals(pfas):  # n:2 Fluorotelomer unsaturated aldehydes
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=C/C=C(\F)C"):
        smiles = smiles.replace("O=C/C=C(\F)", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    elif smiles.startswith("O=CC=C(F)C"):
        smiles = smiles.replace("O=CC=C(F)", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftcas(pfas):  # n:2 Fluorotelomer carboxylic acids
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=C(O)CC"):
        smiles = smiles.replace("O=C(O)C", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftucas(pfas):  # n:2 Fluorotelomer unsaturated carboxylic acids
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=C(O)C=C(F)C"):
        smiles = smiles.replace("O=C(O)C=C(F)", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_ftsas(pfas):  # n:2 Fluorotelomer sulfonic acids
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=S(=O)(O)CCC"):
        smiles = smiles.replace("O=S(=O)(O)CC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_sfas(pfas):  # Semifluorinated n-alkanes
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    atom_num = pfas.get_natom()
    if c_num + f_num + h_num == atom_num and h_num != 0:
        smiles = pfas.get_rdkit_smiles()
        split_smiles = smiles.split("C(", 1)
        if len(split_smiles) == 2:
            alkyl_chain = split_smiles[0]
            pfai_chain = "IC(" + split_smiles[1]
            if determine_alkyl_chain(alkyl_chain) and classify_pfais(PFAS(pfai_chain)):
                return True
    return False


def classify_sfaenes(pfas):  # Semifluorinated n-alkenes
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    atom_num = pfas.get_natom()
    if c_num + f_num + h_num == atom_num and h_num != 0:
        smiles = pfas.get_rdkit_smiles()
        if ("C=CC") in smiles:
            split_smiles = smiles.split("C=C", 1)
            if len(split_smiles) == 2:
                alkyl_chain = split_smiles[0]
                pfai_chain = "I" + split_smiles[1]
                if determine_alkyl_chain(alkyl_chain) and classify_pfais(PFAS(pfai_chain)):
                    return True
        elif ("/C=C/C") in smiles:
            split_smiles = smiles.split("/C=C/", 1)
            if len(split_smiles) == 2:
                alkyl_chain = split_smiles[0]
                pfai_chain = "I" + split_smiles[1]
                if determine_alkyl_chain(alkyl_chain) and classify_pfais(PFAS(pfai_chain)):
                    return True
    return False


def classify_three_acids(pfas):  # n:3 Fluorotelomer carboxylic acids
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=C(O)CCC"):
        smiles = smiles.replace("O=C(O)CC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_unsat_three_acids(pfas):  # n:3 Fluorotelomer unsaturated carboxylic acids
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("O=C(O)C=CC"):
        smiles = smiles.replace("O=C(O)C=C", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    elif smiles.startswith("O=C(O)/C=C/C"):
        smiles = smiles.replace("O=C(O)/C=C/", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classify_n1_ftohs(pfas):  # n:1 Fluorotelomer alcohols
    smiles = pfas.get_rdkit_smiles()
    if smiles.startswith("OCC"):
        smiles = smiles.replace("OC", "I")
        if classify_pfais(PFAS(smiles)):
            return True
    return False


def classifying_fluorotelomers(pfas):
    if classify_ftis(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTIs']
    if classify_ftos(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTOs']
    if classify_ftohs(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTOHs']
    if classify_ftacs(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTACs']
    if classify_ftmacs(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTMACs']
    if classify_monoesters(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 monoPAPs']
    if classify_diesters(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 diPAPs']
    if classify_ftals(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTALs']
    if classify_ftuals(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTUALs']
    if classify_ftcas(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTCAs']
    if classify_ftucas(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTUCAs']
    if classify_ftsas(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:2 FTSAs']
    if classify_sfas(pfas):
        return ['Fluorotelomer PFAA precursors', 'SFAs']
    if classify_sfaenes(pfas):
        return ['Fluorotelomer PFAA precursors', 'SFAenes']
    if classify_three_acids(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:3 Acids']
    if classify_unsat_three_acids(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:3 UAcids']
    if classify_n1_ftohs(pfas):
        return ['Fluorotelomer PFAA precursors', 'n:1 FTOHs']
    return None