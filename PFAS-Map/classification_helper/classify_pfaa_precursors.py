

# PASFs
def classify_pasfs(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    o_num = pfas.get_no()
    s_num = pfas.get_ns()
    atom_num = pfas.get_natom()
    if c_num * 2 + 2 == f_num \
            and o_num == 2 \
            and s_num == 1 \
            and c_num + f_num + s_num + o_num == atom_num \
            and pfas.get_rdkit_smiles().startswith("O=S(=O)(F)"):
        return True
    return False


# PASF-based substances
# Perfluoroalkane sulfonamides (FASAs)
def classify_fasas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    n_num = pfas.get_nn()
    s_num = pfas.get_ns()
    atom_num = pfas.get_natom()
    if c_num * 2 + 1 == f_num \
            and h_num == 2 \
            and o_num == 2 \
            and n_num == 1 \
            and s_num == 1 \
            and c_num + f_num + s_num + o_num + n_num + h_num == atom_num \
            and pfas.get_rdkit_smiles().startswith("NS(=O)(=O)"):
        return True
    return False


# Perfluoroalkanoyl fluorides (PAFs)
def classify_pafs(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if c_num * 2 == f_num \
            and o_num == 1 \
            and c_num + f_num + o_num == atom_num \
            and pfas.get_rdkit_smiles().startswith('O=C(F)'):
        return True
    return False


# Perfluoroalkyl iodides (PFAIs)
def classify_pfais(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    i_num = pfas.get_ni()
    atom_num = pfas.get_natom()
    if c_num * 2 + 1 == f_num \
            and i_num == 1 \
            and c_num + f_num + i_num == atom_num \
            and (pfas.get_rdkit_smiles().endswith("I") or '(I)' in pfas.get_rdkit_smiles()):
        return True
    return False


# Perfluoroalkyl aldehydes and aldehyde hydrates (PFALs)
def classify_pfals(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if (c_num - 1) * 2 + 1 == f_num \
            and h_num == 1 \
            and o_num == 1 \
            and c_num + f_num + o_num + h_num == atom_num \
            and pfas.get_rdkit_smiles().startswith("O=C") \
            and not (pfas.get_rdkit_smiles().startswith("O=C(F)")):
        return True
    return False


def classifying_pfaa_precursors(pfas):
    if classify_pasfs(pfas):
        return ['Perfluoro PFAA precursors', 'PASFs']
    if classify_fasas(pfas):
        return ['Perfluoro PFAA precursors', 'FASAs']
    if classify_pafs(pfas):
        return ['Perfluoro PFAA precursors', 'PAFs']
    if classify_pfais(pfas):
        return ['Perfluoro PFAA precursors', 'PFAIs']
    if classify_pfals(pfas):
        return ['Perfluoro PFAA precursors', 'PFALs']
    return None