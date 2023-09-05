import os
from rdkit import Chem
from rdkit.Chem import AllChem

def worker(protac, warhead):
    protac_mol = Chem.MolFromMolFile(protac)
    warhead_mol = Chem.MolFromMolFile(warhead)

    hit_ats = protac_mol.GetSubstructMatch(warhead_mol)

    warhead_conf = warhead_mol.GetConformer()

    match_dict = {}

    for idx, match_idx in enumerate(hit_ats):
        if idx <= 27:
            match_dict[match_idx] = warhead_conf.GetAtomPosition(idx)


    protac_mol_addh = Chem.AddHs(protac_mol)

    # protac_mol_addh_conf = protac_mol_addh.GetConformer()
    # print(protac_mol_addh_conf.GetAtomPosition(0).x)

    #### first step, This step return 0 is ok, -1 means a bug.
    prbcid = AllChem.EmbedMolecule(protac_mol_addh, coordMap=match_dict)

    #### second step
    atomMap = list(zip(hit_ats, range(len(hit_ats))))
    rmsd = AllChem.AlignMol(protac_mol_addh, warhead_mol, atomMap=atomMap)

    with Chem.SDWriter('./protac_align_to_warhead.sdf') as w:
        w.write(protac_mol_addh)
    print('Done')

    print('Done')




if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-TC/'

    protac = os.path.join(base_dir, 'my_script', 'last_protac_align.sdf')
    warhead = os.path.join(base_dir, 'my_script', 'docked_warhead.sdf')

    worker(protac, warhead)





