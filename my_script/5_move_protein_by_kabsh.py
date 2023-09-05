import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import math
from rdkit.Geometry.rdGeometry import Point3D

def rigid_transform_Kabsch_3D(A, B):
    assert A.shape[1] == B.shape[1]
    num_rows, num_cols = A.shape
    if num_rows != 3:
        raise Exception(f"matrix A is not 3xN, it is {num_rows}x{num_cols}")
    num_rows, num_cols = B.shape
    if num_rows != 3:
        raise Exception(f"matrix B is not 3xN, it is {num_rows}x{num_cols}")


    # find mean column wise: 3 x 1
    centroid_A = np.mean(A, axis=1, keepdims=True)
    centroid_B = np.mean(B, axis=1, keepdims=True)

    # subtract mean
    Am = A - centroid_A
    Bm = B - centroid_B

    H = Am @ Bm.T

    # find rotation
    U, S, Vt = np.linalg.svd(H)

    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        # print("det(R) < R, reflection detected!, correcting for it ...")
        SS = np.diag([1.,1.,-1.])
        R = (Vt.T @ SS) @ U.T
    assert math.fabs(np.linalg.det(R) - 1) < 1e-5

    t = -R @ centroid_A + centroid_B
    return R, t

def get_coord(protac_mol):
    protac_conf = protac_mol.GetConformer()

    coord_list = []
    for idx, atom in enumerate(protac_mol.GetAtoms()):
        coord = np.array(list(protac_conf.GetAtomPosition(idx)))
        coord_list.append(coord)
    return np.stack(coord_list, axis=0)
def worker(protac_src, protac_tgt):
    """
    protac_src: the E3-ligand part is right pose
    protac_tgt: the warhead part is right pose
    """
    protac_src_mol = Chem.MolFromMolFile(protac_src, removeHs=False)
    protac_tgt_mol = Chem.MolFromMolFile(protac_tgt, removeHs=False)

    protac_src_coords = get_coord(protac_src_mol)
    protac_tgt_coords = get_coord(protac_tgt_mol)

    R, t = rigid_transform_Kabsch_3D(protac_src_coords.T, protac_tgt_coords.T)


    test_coord = (R @ protac_src_coords.T + t).T

    protac_src_mol_conf = protac_src_mol.GetConformer()

    for idx, (atom, coord) in enumerate(zip(protac_src_mol.GetAtoms(), test_coord)):
        print(atom, coord)
        point3d = Point3D(*coord)
        protac_src_mol_conf.SetAtomPosition(idx, Point3D(*coord))
        print(protac_src_mol_conf.GetAtomPosition(idx).x)




    with Chem.SDWriter('./last_protac.sdf') as w:
        w.write(protac_src_mol)
    print('Done')

    print('Done')


    # for atom in protac_src_mol.GetAtoms():
    #     print(atom.GetSymbol(), atom.GetIdx())
    # print('---------------------------')
    #
    # for atom in protac_tgt_mol.GetAtoms():
    #     print(atom.GetSymbol(), atom.GetIdx())

if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-TC/'

    protac_src = os.path.join(base_dir, 'my_script', 'last_protac_align.sdf')
    protac_tgt = os.path.join(base_dir, 'my_script', 'protac_align_to_warhead.sdf')



    worker(protac_src, protac_tgt)



