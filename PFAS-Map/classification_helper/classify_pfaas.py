

def classify_pfcas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if (c_num - 1) * 2 + 1 == f_num \
            and h_num == 1 \
            and o_num == 2 \
            and c_num + h_num + f_num + o_num == atom_num \
            and pfas.get_rdkit_smiles().startswith('O=C(O)'):
        return True
    return False


def classify_pfsas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    s_num = pfas.get_ns()
    atom_num = pfas.get_natom()
    if c_num * 2 + 1 == f_num \
            and h_num == 1 \
            and o_num == 3 \
            and s_num == 1 \
            and c_num + h_num + f_num + o_num + s_num == atom_num \
            and pfas.get_rdkit_smiles().startswith('O=S(=O)(O)'):
        return True
    return False


def classify_pfsias(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    s_num = pfas.get_ns()
    atom_num = pfas.get_natom()
    if c_num * 2 + 1 == f_num \
            and h_num == 1 \
            and o_num == 2 \
            and s_num == 1 \
            and c_num + h_num + f_num + o_num + s_num == atom_num \
            and pfas.get_rdkit_smiles().startswith('O=S(O)'):
        return True
    return False


def classify_pfecas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    atom_num = pfas.get_natom()
    if ((c_num - 1) * 2 + 1 == f_num
            and h_num == 1
            and o_num > 2
            and c_num + h_num + f_num + o_num == atom_num
            and pfas.get_rdkit_smiles().startswith('O=C(O)')):
        return True
    return False


def classify_pfesas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    s_num = pfas.get_ns()
    atom_num = pfas.get_natom()
    if (c_num * 2 + 1 == f_num
            and h_num == 1
            and o_num > 3
            and s_num == 1
            and c_num + h_num + f_num + o_num + s_num == atom_num
            and pfas.get_rdkit_smiles().startswith('O=S(=O)(O)')):
        return True
    return False


def classify_pfpas(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    p_num = pfas.get_np()
    atom_num = pfas.get_natom()
    if c_num * 2 + 1 == f_num \
            and h_num == 2 \
            and o_num == 3 \
            and p_num == 1 \
            and c_num + h_num + f_num + o_num + p_num == atom_num \
            and pfas.get_rdkit_smiles().startswith("O=P(O)(O)"):
        return True
    return False


def classify_pfpias(pfas):
    c_num = pfas.get_nc()
    f_num = pfas.get_nf()
    h_num = pfas.get_nh()
    o_num = pfas.get_no()
    p_num = pfas.get_np()
    atom_num = pfas.get_natom()
    if c_num * 2 + 2 == f_num \
            and h_num == 1 \
            and o_num == 2 \
            and p_num == 1 \
            and c_num + h_num + f_num + o_num + p_num == atom_num \
            and pfas.get_rdkit_smiles().startswith("O=P(O)"):
        return True
    return False


def classifying_pfaas(pfas):
    if classify_pfcas(pfas):
        return ['PFAAs', 'PFCAs']
    if classify_pfsas(pfas):
        return ['PFAAs', 'PFSAs']
    if classify_pfsias(pfas):
        return ['PFAAs', 'PFSiAs']
    if classify_pfecas(pfas):
        return ['PFAAs', 'PFECAs']
    if classify_pfesas(pfas):
        return ['PFAAs', 'PFESAs']
    if classify_pfpas(pfas):
        return ['PFAAs', 'PFPAs']
    if classify_pfpias(pfas):
        return ['PFAAs', 'PFPiAs']
    return None

