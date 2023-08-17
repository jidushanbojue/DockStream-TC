from rdkit import Chem
import numpy as np
import os
import math
from rdkit import Chem
from rdkit.Chem import AllChem


# def rigid_transform_Kabsch_3D(A, B):
#     assert A.shape[1] == B.shape[1]
#     num_rows, num_cols = A.shape
#     if num_rows != 3:
#         raise Exception(f"matrix A is not 3xN, it is {num_rows}x{num_cols}")
#     num_rows, num_cols = B.shape
#     if num_rows != 3:
#         raise Exception(f"matrix B is not 3xN, it is {num_rows}x{num_cols}")
#
#     # find mean column wise: 3 x 1
#     centroid_A = np.mean(A, axis=1, keepdims=True)
#     centroid_B = np.mean(B, axis=1, keepdims=True)
#
#     # subtract mean
#     Am = A - centroid_A
#     Bm = B - centroid_B
#
#     H = Am @ Bm.T
#
#     # find rotation
#     U, S, Vt = np.linalg.svd(H)
#
#     R = Vt.T @ U.T
#
#     # special reflection case
#     if np.linalg.det(R) < 0:
#         # print("det(R) < R, reflection detected!, correcting for it ...")
#         SS = np.diag([1.,1.,-1.])
#         R = (Vt.T @ SS) @ U.T
#     assert math.fabs(np.linalg.det(R) - 1) < 1e-5
#
#     t = -R @ centroid_A + centroid_B
#     return R, t
#
# def get_coordinate(docked_pose_sdf, origin_protac_sdf):
#     pass



def align_e3_ligand(docked_E3_ligand_pose, origin_protac):
    docked_mol = Chem.MolFromMolFile(docked_E3_ligand_pose)

    for atom in docked_mol.GetAtoms():
        print(atom.GetIsotope())
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)

    protac_mol = Chem.MolFromMolFile(origin_protac)

    hit_ats = protac_mol.GetSubstructMatch(docked_mol)
    match_dict = {}

    docked_conf = docked_mol.GetConformer()

    for idx, match_idx in enumerate(hit_ats):
        match_dict[match_idx] = docked_conf.GetAtomPosition(idx)

    protac_mol_addh = Chem.AddHs(protac_mol)

    protac_mol_addh_conf = protac_mol_addh.GetConformer()
    print(protac_mol_addh_conf.GetAtomPosition(0).x)

    #### first step
    prbcid = AllChem.EmbedMolecule(protac_mol_addh, coordMap=match_dict, useMacrocycleTorsions=True)

    #### second step
    atomMap = list(zip(hit_ats, range(len(hit_ats))))
    rmsd = AllChem.AlignMol(protac_mol_addh, docked_mol, atomMap=atomMap)

    with Chem.SDWriter('./last_protac_align.sdf') as w:
        w.write(protac_mol_addh)
    print('Done')


def constrained_minimization(docked_E3_ligand_pose, aligned_protac_pose):
    docked_mol = Chem.MolFromMolFile(docked_E3_ligand_pose)

    for atom in docked_mol.GetAtoms():
        print(atom.GetIsotope())
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)

    protac_mol = Chem.MolFromMolFile(aligned_protac_pose)
    hit_ats = protac_mol.GetSubstructMatch(docked_mol)
    print('Done')



if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-TC'
    data_dir = os.path.join(base_dir, 'data/Glide')

    docked_pose_sdf = os.path.join(data_dir, 'docked_E3_ligand.sdf')
    origin_protac_sdf = os.path.join(data_dir, 'crystal_protac_ligprep.sdf')


    aligned_protac_sdf = os.path.join(base_dir, 'my_script', 'last_protac_align.sdf')

    # align_e3_ligand(docked_E3_ligand_pose=docked_pose_sdf,
    #                 origin_protac=origin_protac_sdf)

    constrained_minimization(docked_pose_sdf, aligned_protac_sdf)






