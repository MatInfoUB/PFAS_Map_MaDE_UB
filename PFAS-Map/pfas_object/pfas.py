from rdkit_helper import mol_from_smiles, rdkit_smiles_from_mol
from padel_helper import generate_padel_descriptors, fetch_padel_descriptors
import pandas as pd


class PFAS:

    def __init__(self, input_smiles, exist=False, df_descriptors=pd.DataFrame()):

        self.mol = mol_from_smiles(input_smiles)
        if self.mol is not None:
            self.rdkit_smiles = rdkit_smiles_from_mol(self.mol)
            if not exist:
                try:
                    self.descriptors = generate_padel_descriptors(self.rdkit_smiles)
                except RuntimeError:
                    self.descriptors = {}
            else:
                self.descriptors = fetch_padel_descriptors(self.rdkit_smiles, df_descriptors)
        else:
            self.rdkit_smiles = ""
            self.descriptors = {}

    def get_mol(self):
        return self.mol

    def get_rdkit_smiles(self):
        return self.rdkit_smiles

    def get_descriptors(self):
        return self.descriptors

    def get_nc(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nC"]))
        return -1

    def get_nf(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nF"]))
        return -1

    def get_no(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nO"]))
        return -1

    def get_nh(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nH"]))
        return -1

    def get_nn(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nN"]))
        return -1

    def get_ns(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nS"]))
        return -1

    def get_np(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nP"]))
        return -1

    def get_ni(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nI"]))
        return -1

    def get_natom(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAtom"]))
        return -1

    def get_arom(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAromBond"]))
        return -1