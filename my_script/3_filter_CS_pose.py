import math
import os
from biopandas.pdb import PandasPdb
from rdkit import Chem
import math
import numpy as np

# def distance(p1, p2):
#     x1, y1, z1 = p1
#     x2, y2, z2 = p2
#     return math.hypot(x2-x1, y2-y1, z2-z1)

def distance(p1, p2):
    """
    p1: protac_atom_coord
    p2: protein_atom_coord
    """
    return np.linalg.norm(p1-p2)


def filter_aligned_protacs(protein, aligned_protac_cs, docked_E3_ligand, cutoff):
    """
    protein: E3-liganse
    aligned_protac_cs: Constrained Conformational search result for protacs
    cutoff: 2 angstron
    """

    E3_ligand_mol = Chem.MolFromMolFile(docked_E3_ligand)
    for atom in E3_ligand_mol.GetAtoms():
        atom.SetIsotope(0)

    protein_df = PandasPdb().read_pdb(protein).df['ATOM']
    protein_coord_arr = protein_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    aligned_mols = Chem.SDMolSupplier(aligned_protac_cs)
    hit_ats = aligned_mols[0].GetSubstructMatch(E3_ligand_mol)

    result_protacs = []

    for mol in aligned_mols:
        flag = True
        for idx, atom in enumerate(mol.GetAtoms()):
            if flag == False:
                break
            # xyz = list(mol.GetConformer().GetAtomPosition(1))

            if idx in hit_ats:
                continue

            protac_atom_coord = np.array(mol.GetConformer().GetAtomPosition(idx))

            for protein_atom_coord in protein_coord_arr:
                dist = distance(protac_atom_coord, protein_atom_coord)

                if dist < cutoff:
                    flag = False
                    break
        if flag:
            result_protacs.append(mol)
    print(len(result_protacs))
    print('Done')



if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-TC/'

    protein = os.path.join(base_dir, 'data', 'CSearch_2', 'E3_ligase.pdb')
    aligned_protac_cs = os.path.join(base_dir, 'data', 'CSearch_2', 'aligned_protac_cs.sdf')
    docked_E3_ligand = os.path.join(base_dir, 'data', 'CSearch_2', 'docked_E3_ligand.sdf')

    filter_aligned_protacs(protein, aligned_protac_cs, docked_E3_ligand, 4.1)




