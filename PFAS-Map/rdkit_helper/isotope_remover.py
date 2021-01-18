from rdkit import Chem
"""Source: https://sourceforge.net/p/rdkit/mailman/message/36877847/"""


def MolWithoutIsotopesToSmiles(mol):
   atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
   for atom, isotope in atom_data:
       if isotope:
           atom.SetIsotope(0)
   return mol


