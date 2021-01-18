

# Perfluorinated alkanes (PFAs)
def classify_pfas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    atom_num = pfas.get_natom()
    if atom_num == c_num + f_num and c_num * 2 + 2 == f_num:
        return True
    return False


# Perfluorinated alkenes (PFAenes)
def classify_pfaenes(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    atom_num = pfas.get_natom()
    if atom_num == c_num + f_num and \
            c_num * 2 == f_num \
            and '=' in pfas.get_rdkit_smiles():
        return True
    return False


# Perfluoroalkyl alcohols (PFACs)
def classify_pfacs(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if atom_num == c_num + f_num + o_num + h_num \
            and h_num == 1 \
            and o_num == 1 \
            and pfas.get_rdkit_smiles().startswith('OC'):
        return True
    return False


# Perfluoroalkyl ketones (PFAKs)
def classify_pfaks(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if c_num * 2 == f_num and atom_num == c_num + f_num + o_num \
            and o_num == 1 \
            and pfas.get_rdkit_smiles().startswith('O=C'):
        return True
    return False


def classifying_other_perfluoroalkyls(pfas):
    if classify_pfas(pfas):
        return ['Non-PFAA perfluoroalkyls', 'PFAs']
    if classify_pfaenes(pfas):
        return ['Non-PFAA perfluoroalkyls', 'PFAenes']
    if classify_pfacs(pfas):
        return ['Non-PFAA perfluoroalkyls', 'PFACs']
    if classify_pfaks(pfas):
        return ['Non-PFAA perfluoroalkyls', 'PFAKs']
    return None