import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem

def connect_molecules(smiles1, smiles2, output_name, virtual_symbol="*", nbonds=1):
    """
    连接两个带虚粒子的分子，通过虚粒子的邻接原子生成整体的分子文件。
    Args:
    - smiles1 (str): 第一个带虚粒子的 SMILES。
    - smiles2 (str): 第二个带虚粒子的 SMILES。
    - output_name (str): 输出文件的前缀名称。
    - virtual_symbol (str): 虚粒子的符号。

    Returns:
    - str: 输出的整体 SMILES。
    """
    # 将 SMILES 转换为 RDKit 分子对象
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    nbonds = int(nbonds)
    
    if not mol1 or not mol2:
        raise ValueError("Failed to parse one or both SMILES strings. Check input.")
    
    # 找到第一个分子中虚粒子的邻接原子
    atom1 = None
    for atom in mol1.GetAtoms():
        if atom.GetSymbol() == virtual_symbol:
            neighbors = atom.GetNeighbors()
            if neighbors:
                atom1 = neighbors[0].GetIdx()
            break
    
    # 找到第二个分子中虚粒子的邻接原子
    atom2 = None
    for atom in mol2.GetAtoms():
        if atom.GetSymbol() == virtual_symbol:
            neighbors = atom.GetNeighbors()
            if neighbors:
                atom2 = neighbors[0].GetIdx()
            break
    
    if atom1 is None or atom2 is None:
        raise ValueError("Failed to find neighboring atoms for virtual atoms in one or both molecules.")
    
    # 合并分子
    combined = Chem.CombineMols(mol1, mol2)
    rw_mol = Chem.RWMol(combined)
    
    # 调整索引偏移
    offset = mol1.GetNumAtoms()  # 第二个分子的原子索引需要偏移
    atom2 = atom2 + offset

    if nbonds == 1:
        # 添加键
        rw_mol.AddBond(atom1, atom2, Chem.BondType.SINGLE)
    elif nbonds == 2:
        # 添加键
        rw_mol.AddBond(atom1, atom2, Chem.BondType.DOUBLE)
    elif nbonds == 3:
        # 添加键
        rw_mol.AddBond(atom1, atom2, Chem.BondType.TRIPLE)
    
    # 移除虚粒子
    for atom in mol1.GetAtoms():
        if atom.GetSymbol() == virtual_symbol:
            rw_mol.RemoveAtom(atom.GetIdx())
    for atom in mol2.GetAtoms():
        if atom.GetSymbol() == virtual_symbol:
            rw_mol.RemoveAtom(atom.GetIdx() + offset - 1)  # Adjust for removed atoms

    # 更新分子并重新构建
    final_mol = Chem.Mol(rw_mol)
    Chem.SanitizeMol(final_mol)  # 确保分子结构完整
    
    # 转换回 SMILES
    combined_smiles = Chem.MolToSmiles(final_mol)
    with open(f"{output_name}_connected.smi", "w") as f:
        f.write(combined_smiles)
    
    return combined_smiles

def process_smiles(smiles, idx=1):
    """
    Removes the dummy atom [*] from the SMILES string and adds an atom mapping number to the connected atom.
    Returns the processed SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"Failed to parse SMILES: {smiles}")
        return None

    # Find the dummy atom [*]
    dummy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if not dummy_atoms:
        print("No dummy atom [*] found in the SMILES")
        return smiles

    if len(dummy_atoms) > 1:
        print("Multiple dummy atoms found, only processing the first one")

    dummy_atom = dummy_atoms[0]
    dummy_idx = dummy_atom.GetIdx()
    # Get the neighbor atom of the dummy atom
    neighbors = dummy_atom.GetNeighbors()
    if not neighbors:
        print("Dummy atom has no neighbors")
        return None

    # Assuming we want to add the mapping to the connected atom
    connected_atom = neighbors[0]

    # **Set the atom mapping number before modifying the molecule**
    connected_atom.SetAtomMapNum(idx)

    # Create an editable molecule
    mol_edit = Chem.RWMol(mol)

    # **Remove the bond between the dummy atom and the connected atom**
    mol_edit.RemoveBond(dummy_idx, connected_atom.GetIdx())

    # **Remove the dummy atom**
    mol_edit.RemoveAtom(dummy_idx)

    # Get the processed molecule
    processed_mol = mol_edit.GetMol()

    # **Generate the SMILES without canonicalization to preserve atom order**
    processed_smiles = Chem.MolToSmiles(processed_mol, isomericSmiles=True, canonical=False)

    return processed_smiles

def read_molecule_from_pdb(pdb, template_smiles=None):
    """从PDB格式的文件中读取分子，并尝试重建键信息"""
    mol = Chem.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
    if template_smiles:
        ref_mol = Chem.MolFromSmiles(template_smiles)
        ref_mol = Chem.AddHs(ref_mol)
        AllChem.AssignBondOrdersFromTemplate(ref_mol, mol)
    if mol:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                         Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                         Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                         Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
    return mol

def get_mapping_name_from_pdb(pdb, template_smiles=None):
    """读取 PDB 文件中的分子，并尝试重建键序信息。返回 mol 对象和映射编号对应的 PDB 原子名称字典。"""
    # 读取 PDB 文件中的分子
    mol = Chem.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
    if not mol:
        print(f"Failed to load molecule from PDB file: {pdb}")
        return None, None

    mapping_atom_names = {}  # 用于存储映射编号对应的 PDB 原子名称

    if template_smiles:
        # 生成带有映射编号的参考分子
        ref_mol = Chem.MolFromSmiles(template_smiles)
        ref_mol = Chem.AddHs(ref_mol)

        # 使用参考分子分配键序
        mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol)

        # 转移原子映射编号，并获取对应的 PDB 原子名称
        for ref_atom, mol_atom in zip(ref_mol.GetAtoms(), mol.GetAtoms()):
            if ref_atom.HasProp("molAtomMapNumber"):
                map_num = int(ref_atom.GetProp("molAtomMapNumber"))
                mol_atom.SetAtomMapNum(map_num)
                pdb_info = mol_atom.GetPDBResidueInfo()
                if pdb_info:
                    atom_name = pdb_info.GetName().strip()
                    mapping_atom_names[map_num] = atom_name
            else:
                mol_atom.SetAtomMapNum(0)

    if mol:
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        Chem.SanitizeMol(mol)

    return mapping_atom_names

def find_smiles_pattern(pdb, smiles="NCC(O)O", return_first=False):
    """Finds SMILES patterns in a PDB file."""
    mol = read_molecule_from_pdb(pdb)
    if not mol:
        return "Failed to load molecule from PDB"

    pattern = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(pattern)  # Ensure the pattern is sanitized for matching
    matches = mol.GetSubstructMatches(pattern)

    atom_names = []
    if matches:
        # Output all matching atom names
        for match in matches:
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            atom_names.extend([atom.GetPDBResidueInfo().GetName() for atom in atoms])
            if return_first:
                break
        return atom_names
    else:
        return []

def map_to_amber_mc(molecule, 
                    matched_atoms,
                    head_name=None,
                    tail_name=None,
                    main_chain=None,
                    omit_name=None,
                    pre_head_type='C',
                    pre_tail_type='N',
                    charge='0.0'):
    """
    将匹配的原子基于它们的化学环境映射到Amber mc文件的HEAD_NAME、TAIL_NAME和MAIN_CHAIN字段。
    Args:
    molecule (rdkit.Chem.Mol): RDKit分子对象。
    matched_atoms (list of tuples): 每个元组包含匹配的原子名称。
    head_name (str): 头部原子的名称。
    tail_name (str): 尾部原子的名称。
    main_chain (str): 主链原子的名称。
    omit_name (list of str): 需要忽略的原子名称。
    pre_head_type (str): 头部原子的类型。
    pre_tail_type (str): 尾部原子的类型。
    charge (str): 分子的电荷。

    Returns:
    dict: 包含HEAD_NAME, TAIL_NAME和MAIN_CHAIN的字典。
    """

    amber_mapping = {
        'HEAD_NAME': [] if head_name is None else head_name,
        'TAIL_NAME': [] if tail_name is None else tail_name,
        'MAIN_CHAIN': [] if main_chain is None else main_chain,
        'OMIT_NAME': [] if omit_name is None else omit_name,
        'PRE_HEAD_TYPE': pre_head_type,
        'PRE_TAIL_TYPE': pre_tail_type,
        'CHARGE': charge
    }

    O_mappings = []
    for atom in molecule.GetAtoms():
        if atom.GetPDBResidueInfo().GetName() not in matched_atoms:
            continue

        # N原子的度为3，相邻原子为C,和2个H
        if atom.GetSymbol() == 'N' and atom.GetTotalDegree() == 3 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']) == 2 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']) == 1:
            # amber_mapping['HEAD_NAME'] = atom.GetPDBResidueInfo().GetName().strip()
            amber_mapping['HEAD_NAME'].append(atom.GetPDBResidueInfo().GetName().strip())

        # C原子的度为4, 相邻原子为=O，-O，C
        elif atom.GetSymbol() == 'C' and \
            atom.GetTotalDegree() == 4 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) == 2 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']) == 1:
            # amber_mapping['TAIL_NAME'] = atom.GetPDBResidueInfo().GetName().strip()
            amber_mapping['TAIL_NAME'].append(atom.GetPDBResidueInfo().GetName().strip())

        # 链接C和N的原子为主链原子
        if atom.GetSymbol() == 'C' and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'N']) == 1 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) == 0:
            # amber_mapping['MAIN_CHAIN'] = atom.GetPDBResidueInfo().GetName().strip()
            amber_mapping['MAIN_CHAIN'].append(atom.GetPDBResidueInfo().GetName().strip())

    # 开始获取OMIT_NAME, 
    # 1. 与TAIL_NAME相邻的OH原子，包括O原子及其相邻的H原子
    for atom in molecule.GetAtoms():
        if atom.GetPDBResidueInfo().GetName().strip() in amber_mapping['TAIL_NAME']:
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 1:
                    amber_mapping['OMIT_NAME'].append(n.GetPDBResidueInfo().GetName().strip())
                    amber_mapping['OMIT_NAME'].extend([n.GetPDBResidueInfo().GetName().strip() for n in n.GetNeighbors() if n.GetSymbol() == 'H'])

                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 0:
                    O_mapping = n.GetPDBResidueInfo().GetName().strip() # O_mapping是O原子的名称
                    O_mappings.append(O_mapping)

    for atom in molecule.GetAtoms():
        if atom.GetPDBResidueInfo().GetName().strip() == amber_mapping['HEAD_NAME'][0]:
            # 2. 与HEAD_NAME相邻的NH原子，有2个H，但是只删除一个H
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'H' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'N']) == 1:
                    amber_mapping['OMIT_NAME'].append(n.GetPDBResidueInfo().GetName().strip())
                    break

    return amber_mapping, O_mappings

def modify_ac(ac, atom_map, O_mappings):
    '''
    Modifies the ac file using sed commands.
    '''
    # 1. Modify N3 to N
    head_names = atom_map['HEAD_NAME']
    with open(ac, 'r') as f:
        lines = f.readlines()
    for head_name in head_names:
        for i, line in enumerate(lines):
            # 列14开始
            len_head_name = len(head_name)
            if head_name in line[13:13+len_head_name]:
                lines[i] = line.replace(head_name, 'N'.ljust(len_head_name))

    # 2. Modify C6 to C
    tail_names = atom_map['TAIL_NAME']
    for tail_name in tail_names:
        for i, line in enumerate(lines):
            len_tail_name = len(tail_name)
            if tail_name in line[13:13+len_tail_name]:
                lines[i] = line.replace(tail_name, 'C'.ljust(len_tail_name))

    for O_mapping in O_mappings:
        # 3. Modify O1 to O
        for i, line in enumerate(lines):
            len_O_mapping = len(O_mapping)
            if O_mapping in line[13:13+len_O_mapping]:
                lines[i] = line.replace(O_mapping, 'O'.ljust(len_O_mapping))

    # 4. Modify C5 to CA
    main_chains = atom_map['MAIN_CHAIN']
    for main_chain in main_chains:
        for i, line in enumerate(lines):
            len_main_chain = len(main_chain)
            if main_chain in line[13:13+len_main_chain]:
                lines[i] = line.replace(main_chain, 'CA'.ljust(len_main_chain))

    with open(ac, 'w') as f:
        f.writelines(lines)

def modify_prepin(prepin, atom_map, O_mappings):
    '''
    Modifies the prepin file using sed commands.
    '''
    # 1. Modify N3 to N
    head_names = atom_map['HEAD_NAME']
    for head_name in head_names:
        if len(head_name) == 2:
            head_name = f'{head_name} '
        cmd = f"sed -i 's/{head_name}   NT/N     N /g' {prepin}"
        print(cmd)
        os.system(cmd)

    # 2. Modify C6 to C
    tail_names = atom_map['TAIL_NAME']
    for tail_name in tail_names:
        if len(tail_name) == 2:
            tail_name = f'{tail_name} '
        cmd = f"sed -i 's/{tail_name}   C/C     C/g' {prepin}"
        print(cmd)
        os.system(cmd)

    for O_mapping in O_mappings:
        # 3. Modify O1 to O
        if len(O_mapping) == 2:
            O_mapping = f'{O_mapping} '
        elif len(O_mapping) == 3:
            O_mapping = f'{O_mapping}'
        elif len(O_mapping) == 1:
            O_mapping = f'{O_mapping}   '
        cmd = f"sed -i 's/{O_mapping}   O/O     O/g' {prepin}"
        print(cmd)
        os.system(cmd)

    # 4. Modify C5 to CA
    main_chains = atom_map['MAIN_CHAIN']
    for main_chain in main_chains:
        if len(main_chain) == 2:
            main_chain = f'{main_chain} '
        cmd = f"sed -i 's/{main_chain}   CT/CA    CT/g' {prepin}"
        print(cmd)
        os.system(cmd)

def write_amber_mc(amber_map, output_file):
    """
    Writes the Amber mapping to an mc file.
    """
    with open(f'{output_file}.mc', 'w') as f:
        f.write(f'HEAD_NAME {amber_map["HEAD_NAME"][0]}\n')
        f.write(f'TAIL_NAME {amber_map["TAIL_NAME"][0]}\n')
        f.write(f'MAIN_CHAIN {amber_map["MAIN_CHAIN"][0]}\n')
        for omit_name in amber_map["OMIT_NAME"]:
            f.write(f'OMIT_NAME {omit_name}\n')
        f.write(f'PRE_HEAD_TYPE {amber_map["PRE_HEAD_TYPE"]}\n')
        f.write(f'PRE_TAIL_TYPE {amber_map["PRE_TAIL_TYPE"]}\n')
        f.write(f'CHARGE {amber_map["CHARGE"]}\n')
    return f'{output_file}.mc'

def smi_to_sdf(smi, name=None):
    """Converts SMILES to 3D coordinates and handles atom mapping numbers."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print(f"Failed to parse SMILES: {smi}")
        return None

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Set atom names based on atom mapping numbers
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            atom_name = f"Map{map_num}"
            # Ensure atom name is at most 4 characters for PDB format
            atom_name = atom_name[:4]
        else:
            atom_name = atom.GetSymbol()
        pdb_info = Chem.AtomPDBResidueInfo()
        pdb_info.SetName(atom_name.ljust(4))
        atom.SetPDBResidueInfo(pdb_info)

    # Write to SDF
    writer = Chem.SDWriter(f'{name}.sdf')
    writer.write(mol)
    writer.close()

    return f'{name}.sdf'

def sdf_to_ac(sdf, name=None):
    """将SDF转换为Amber mc文件"""
    '''
    antechamber -fi sdf -i R1A.sdf -bk R1A -fo ac -o R1A.ac -c bcc -at amber
    '''
    cmd = f'antechamber -fi sdf -i {sdf} -bk {name} -fo ac -o {name}.ac -c bcc -at amber'
    print(cmd)
    os.system(cmd)
    return f'{name}.ac'

def ac_to_pdb(ac, name=None):
    """将Amber mc文件转换为PDB文件"""
    '''
    antechamber -fi ac -i R1A.ac -bk R1A -fo pdb -o R1A.pdb
    '''
    cmd = f'antechamber -fi ac -i {ac} -bk {name} -fo pdb -o {name}.pdb'
    print(cmd)
    os.system(cmd)
    return f'{name}.pdb'

def mc_to_prepin(ac, mc, prepin=None, name=None):
    """Converts Amber MC file to prepin file."""
    if not name:
        name = mc.replace('.mc', '')
    cmd = f'prepgen -i {ac} -o {prepin} -m {mc} -rn {name}'
    print(cmd)
    os.system(cmd)
    return prepin

def prep_to_frcmod(prep, frcmod=None):
    """Converts prepin file to frcmod file."""
    cmd = f'parmchk2 -i {prep} -f prepi -o {frcmod} -a Y'
    print(cmd)
    os.system(cmd)

def ac_to_frcmod(ac, frcmod=None):
    """Converts ac file to frcmod file."""
    cmd = f'parmchk2 -i {ac} -f ac -o {frcmod} -a Y'
    print(cmd)
    os.system(cmd)

def get_omit_names(mol, placeholder_atoms, nbonds=1):
    omit_names = []
    # Include hydrogens attached to the placeholder atoms (if any)
    i = 0
    for atom in mol.GetAtoms():
        print('atom.GetSymbol():', atom.GetNeighbors()[0].GetPDBResidueInfo().GetName().strip())
        # 是H，且它的邻原子是placeholder_atoms
        if atom.GetSymbol() == 'H' and atom.GetNeighbors()[0].GetPDBResidueInfo().GetName().strip() == placeholder_atoms:
            omit_name = atom.GetPDBResidueInfo().GetName().strip()
            omit_names.append(omit_name)
            i += 1
            if i == int(nbonds):
                break
    return omit_names

def merge_frcmod(frcmod_list, output_file):
    """Merges multiple frcmod files into a single file."""
    pass

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Generate non-standard amino acid to Amber mc file')
    parser.add_argument('--smiles1', help='First SMILES string of the non-standard amino acid', required=True)
    parser.add_argument('--smiles2', help='Second SMILES string of the non-standard amino acid', required=True)
    parser.add_argument('--nbonds', help='Number of bonds to break between the two SMILES strings', default=1)
    parser.add_argument('--name1', help='First name of the non-standard amino acid', default='R1A')
    parser.add_argument('--name2', help='Second name of the non-standard amino acid', default='R2A')
    parser.add_argument('--charge1', help='First charge of the non-standard amino acid', default='0')
    parser.add_argument('--charge2', help='Second charge of the non-standard amino acid', default='0')
    parser.add_argument('--output', help='Output directory', default=None)
    args = parser.parse_args()

    print(f'args.smiles1: {args.smiles1}')
    print(f'args.smiles2: {args.smiles2}')

    # Ensure the name is 3 characters
    if len(args.name1) != 3 or (args.smiles2 and len(args.name2) != 3):
        raise ValueError('Name of the non-standard amino acid must be 3 characters')
    
    args.name1 = args.name1.upper()

    if not args.output:
        args.output = args.name1

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        raise ValueError(f'{args.output} already exists')
    
    # Change to output directory
    pwd = os.getcwd()
    os.chdir(args.output)

    # Generate the SDF file
    pre_smiles1 = process_smiles(args.smiles1)
    pre_smiles2 = process_smiles(args.smiles2)
    sdf1 = smi_to_sdf(pre_smiles1, args.name1)
    sdf2 = smi_to_sdf(pre_smiles2, args.name2)
    if sdf1 is None or sdf2 is None:
        sys.exit(1)

    # Generate the AC file
    ac1 = sdf_to_ac(sdf1, args.name1)
    ac2 = sdf_to_ac(sdf2, args.name2)

    # Generate the PDB file
    pdb1 = ac_to_pdb(ac1, args.name1)
    pdb2 = ac_to_pdb(ac2, args.name2)
    mol1 = read_molecule_from_pdb(pdb1, template_smiles=pre_smiles1)
    mol2 = read_molecule_from_pdb(pdb2, template_smiles=pre_smiles2)

    # Read the molecule from the PDB file
    placeholder_atom_names_1 = get_mapping_name_from_pdb(pdb1, template_smiles=pre_smiles1)
    placeholder_atom_names_2 = get_mapping_name_from_pdb(pdb2, template_smiles=pre_smiles2)
    placeholder_atoms_1 = placeholder_atom_names_1[1]
    placeholder_atoms_2 = placeholder_atom_names_2[1]

    # write the covalent atom names
    with open('covalent_info.csv', 'w') as f:
        f.write('Res1,Atom1,Res2,Atom2,nBonds,Frcmod\n')
        for _, name1 in placeholder_atom_names_1.items():
            for _, name2 in placeholder_atom_names_2.items():
                frcmod = f'frcmod.{args.name1}_{args.name2}'
                frcmod = os.path.join(pwd, args.output, frcmod)
                f.write(f'{args.name1},{name1},{args.name2},{name2},{args.nbonds},{frcmod}\n')
                break
    
    # Proceed with your existing functions, using placeholder_atom_names where necessary
    matched_atoms_1 = find_smiles_pattern(pdb1, return_first=True)
    matched_atoms_2 = find_smiles_pattern(pdb2, return_first=True)
    
    omit_names_1 = get_omit_names(mol1, placeholder_atoms_1, nbonds=args.nbonds)
    omit_names_2 = get_omit_names(mol2, placeholder_atoms_2, nbonds=args.nbonds)

    # Map to Amber mc
    amber_map_1, O_mapping_1 = map_to_amber_mc(mol1, matched_atoms_1, charge=args.charge1, omit_name=omit_names_1)
    amber_map_2, O_mapping_2 = map_to_amber_mc(mol2, matched_atoms_2, charge=args.charge1, omit_name=omit_names_2)

    mc1 = write_amber_mc(amber_map_1, args.name1)
    mc2 = write_amber_mc(amber_map_2, args.name2)
    prepin1 = mc_to_prepin(ac1, mc1, f'{args.name1}.prepin')
    prepin2 = mc_to_prepin(ac2, mc2, f'{args.name2}.prepin')
    modify_prepin(prepin1, amber_map_1, O_mapping_1)
    modify_prepin(prepin2, amber_map_2, O_mapping_2)
    prep_to_frcmod(prepin1, f'frcmod.{args.name1}')
    prep_to_frcmod(prepin2, f'frcmod.{args.name2}')

    # combine the smi1 and smi2
    combined_smiles = connect_molecules(args.smiles1, args.smiles2, f'{args.name1}_{args.name2}', nbonds=args.nbonds) # 生成连接后的SMILES
    combined_sdf = smi_to_sdf(process_smiles(combined_smiles), f'{args.name1}_{args.name2}') # 生成连接后的SDF
    combined_ac = sdf_to_ac(combined_sdf, f'{args.name1}_{args.name2}') # 生成连接后的AC
    combined_pdb = ac_to_pdb(combined_ac, f'{args.name1}_{args.name2}') # 生成连接后的PDB
    combined_mol = read_molecule_from_pdb(combined_pdb)
    combined_matched_atoms = find_smiles_pattern(combined_pdb)
    combined_charges = float(args.charge1) + float(args.charge2)
    combined_amber_map, combined_O_mapping = map_to_amber_mc(combined_mol, combined_matched_atoms, charge=combined_charges)
    modify_ac(combined_ac, combined_amber_map, combined_O_mapping)
    ac_to_frcmod(combined_ac, f'frcmod.{args.name1}_{args.name2}')

    # Change back to the original directory
    os.chdir(pwd)

    ## Example usage:
    # python generate_covalent_template.py --smiles1 "[*]C(C)C[C@@H](N)C(=O)O" --name1 R1A --smiles2 "[*]CC(C)C[C@@H](N)C(=O)O" --name2 R2A --output R1A_R2A --nbonds 2