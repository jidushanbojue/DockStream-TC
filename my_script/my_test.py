from rdkit import Chem
import os
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AllChem import AlignMol, GetAlignmentTransform

# smi = 'c1cc(F)ccc1Cl'
# mol = Chem.MolFromSmiles(smi)
#
# patt = Chem.MolFromSmarts('ClccccF')
#
# hit_ats = list(mol.GetSubstructMatch(patt))
# print('Done')
def align_e3_ligand(docked_E3_ligand_pose, origin_protac):
    docked_mol = Chem.MolFromMolFile(docked_E3_ligand_pose)

    for atom in docked_mol.GetAtoms():
        print(atom.GetIsotope())
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)


    protac_mol = Chem.MolFromMolFile(origin_protac)

    print(Chem.MolToSmiles(docked_mol))
    print(Chem.MolToSmiles(protac_mol))

    # for idx, atom in enumerate(docked_mol.GetAtoms()):
    #     atom.SetAtomMapNum(idx)
    # print(Chem.MolToSmiles(docked_mol))
    #
    # with Chem.SDWriter('protac_mol_atommap.sdf') as w:
    #     w.write(protac_mol)

    # d = rdMolDraw2D.MolDraw2DCairo(400, 400)
    # d.DrawMolecule(docked_mol)
    # d.WriteDrawingText('./atom_annotation_1.png')


    # for idx, atom in enumerate(protac_mol.GetAtoms()):
    #     print(atom.GetIdx())
    #     atom.SetAtomMapNum(idx)
    # print(Chem.MolToSmiles(protac_mol))

    hit_ats = protac_mol.GetSubstructMatch(docked_mol)

    match_dict = {}

    # print(len(docked_mol.GetAtoms()))

    docked_conf = docked_mol.GetConformer()

    for idx, match_idx in enumerate(hit_ats):
        match_dict[match_idx] = docked_conf.GetAtomPosition(idx)


    protac_mol_addh = Chem.AddHs(protac_mol)

    protac_mol_addh_conf = protac_mol_addh.GetConformer()
    print(protac_mol_addh_conf.GetAtomPosition(0).x)

    prbcid = AllChem.EmbedMolecule(protac_mol_addh, coordMap=match_dict, useMacrocycleTorsions=True)


    print(Chem.MolToMolBlock(protac_mol_addh))
    print(len(protac_mol_addh.GetAtoms()))

    with Chem.SDWriter('./result_align.sdf') as w:
        w.write(protac_mol_addh)


    # protac_conf = protac_mol_addh.GetConformer()
    #
    # print(protac_conf.GetAtomPosition(0))



    print('Done')



if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-TC'
    data_dir = os.path.join(base_dir, 'data/Glide')

    docked_pose_sdf = os.path.join(data_dir, 'docked_E3_ligand.sdf')
    origin_protac_sdf = os.path.join(data_dir, 'crystal_protac_ligprep.sdf')

    # result_sdf = os.path.join(data_dir, 'result_align.sdf')

    result_sdf = origin_protac_sdf

    # align_e3_ligand(docked_E3_ligand_pose=docked_pose_sdf, origin_protac=origin_protac_sdf)

    docked_mol = Chem.MolFromMolFile(docked_pose_sdf)
    for atom in docked_mol.GetAtoms():
        print(atom.GetIsotope())
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)

    result_mol = Chem.MolFromMolFile(result_sdf)
    # rms = AlignMol(docked_mol, result_mol)

    hit_ats = result_mol.GetSubstructMatch(docked_mol)

    atomMap = list(zip(hit_ats, range(len(hit_ats))))
    # atomMap = list(zip(range(len(hit_ats)), hit_ats))

    # result = GetAlignmentTransform(result_mol, docked_mol, atomMap=atomMap)

    result = AlignMol(result_mol, docked_mol, atomMap=atomMap)


    print(Chem.MolToMolBlock(result_mol))
    print(Chem.MolToMolBlock(docked_mol))

    with Chem.SDWriter('./result_mol_crystal_protac_ligprep.sdf') as w:
        w.write(result_mol)




    print('Done')










