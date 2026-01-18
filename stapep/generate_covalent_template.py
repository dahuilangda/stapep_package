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

def _smiles_has_atom_map(smiles):
    if not smiles:
        return False
    return re.search(r':\d+\]', smiles) is not None

def _strip_atom_map(smiles):
    if not smiles:
        return smiles
    return re.sub(r':\d+\]', ']', smiles)

def _relax_bracket_hcount(smiles):
    if not smiles:
        return smiles
    def repl(match):
        content = match.group(1)
        content = re.sub(r'H(?![a-z])\d*', '', content)
        return f'[{content}]'
    return re.sub(r'\[([^\]]+)\]', repl, smiles)

def _atom_name(atom):
    info = atom.GetPDBResidueInfo()
    if info:
        return info.GetName().strip()
    return atom.GetSymbol()

def _is_head_like_n(atom):
    if atom.GetSymbol() != 'N':
        return False
    heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
    c_neighbors = [n for n in heavy_neighbors if n.GetSymbol() == 'C']
    return len(heavy_neighbors) == 1 and len(c_neighbors) == 1

def _is_carboxyl_carbon(atom):
    return atom.GetSymbol() == 'C' and len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) == 2

def _find_ca_tail_from_head(head_atom, in_scope):
    for c in head_atom.GetNeighbors():
        if c.GetSymbol() != 'C' or not in_scope(c):
            continue
        if _is_carboxyl_carbon(c):
            continue
        for nb in c.GetNeighbors():
            if not in_scope(nb):
                continue
            if _is_carboxyl_carbon(nb):
                return c, nb
    return None, None

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
    template_smiles = _strip_atom_map(smiles) if _smiles_has_atom_map(smiles) else None
    try:
        mol = read_molecule_from_pdb(pdb, template_smiles=template_smiles)
    except Exception:
        mol = read_molecule_from_pdb(pdb)
    if not mol:
        return "Failed to load molecule from PDB"

    pattern_smiles = template_smiles or smiles
    pattern = Chem.MolFromSmiles(pattern_smiles)
    if not pattern:
        return []
    Chem.SanitizeMol(pattern)  # Ensure the pattern is sanitized for matching
    matches = mol.GetSubstructMatches(pattern)
    if not matches and _smiles_has_atom_map(smiles):
        relaxed_smiles = _relax_bracket_hcount(pattern_smiles)
        if relaxed_smiles != pattern_smiles:
            relaxed_pattern = Chem.MolFromSmiles(relaxed_smiles)
            if relaxed_pattern:
                Chem.SanitizeMol(relaxed_pattern)
                matches = mol.GetSubstructMatches(relaxed_pattern)

    atom_names = []
    if matches:
        pattern_with_map = Chem.MolFromSmiles(smiles)
        mapped_atom_indices = []
        if pattern_with_map:
            mapped_atom_indices = [
                idx for idx, atom in enumerate(pattern_with_map.GetAtoms())
                if atom.GetAtomMapNum() > 0
            ]
        use_atom_map = len(mapped_atom_indices) > 0
        # Output all matching atom names
        for match in matches:
            if use_atom_map:
                atoms = [mol.GetAtomWithIdx(match[idx]) for idx in mapped_atom_indices]
                atom_names.extend([atom.GetPDBResidueInfo().GetName().strip() for atom in atoms])
                break
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            atom_names.extend([atom.GetPDBResidueInfo().GetName().strip() for atom in atoms])
            if return_first:
                break
        return list(dict.fromkeys(atom_names))
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
    matched_atom_names = {name.strip() for name in matched_atoms} if matched_atoms else set()
    all_atom_names = {_atom_name(atom) for atom in molecule.GetAtoms()}
    use_match_filter = bool(matched_atom_names)
    if use_match_filter and len(matched_atom_names) >= len(all_atom_names):
        use_match_filter = False

    def _resolve_mapping(use_filter):
        def in_scope(atom):
            if not use_filter:
                return True
            return _atom_name(atom) in matched_atom_names

        head_atom = None
        main_chain_atom = None
        tail_atom = None

        head_candidates = [
            atom for atom in molecule.GetAtoms()
            if in_scope(atom) and _is_head_like_n(atom)
        ]

        for head in head_candidates:
            ca, tail = _find_ca_tail_from_head(head, in_scope)
            if ca and tail:
                head_atom = head
                main_chain_atom = ca
                tail_atom = tail
                break

        if head_atom is None and head_candidates:
            head_atom = head_candidates[0]

        if head_atom and (main_chain_atom is None or tail_atom is None):
            ca, tail = _find_ca_tail_from_head(head_atom, in_scope)
            if main_chain_atom is None:
                main_chain_atom = ca
            if tail_atom is None:
                tail_atom = tail

        if head_atom is None or main_chain_atom is None or tail_atom is None:
            for atom in molecule.GetAtoms():
                if not in_scope(atom) or not _is_carboxyl_carbon(atom):
                    continue
                for ca in atom.GetNeighbors():
                    if ca.GetSymbol() != 'C' or not in_scope(ca):
                        continue
                    if _is_carboxyl_carbon(ca):
                        continue
                    n_neighbors = [n for n in ca.GetNeighbors() if in_scope(n) and _is_head_like_n(n)]
                    if n_neighbors:
                        if head_atom is None:
                            head_atom = n_neighbors[0]
                        if main_chain_atom is None:
                            main_chain_atom = ca
                        if tail_atom is None:
                            tail_atom = atom
                        break
                if head_atom and main_chain_atom and tail_atom:
                    break

        if main_chain_atom is None and head_atom and tail_atom:
            for atom in molecule.GetAtoms():
                if not in_scope(atom):
                    continue
                if atom.GetSymbol() != 'C':
                    continue
                if head_atom in atom.GetNeighbors() and tail_atom in atom.GetNeighbors():
                    main_chain_atom = atom
                    break

        if tail_atom is None and main_chain_atom:
            for nb in main_chain_atom.GetNeighbors():
                if in_scope(nb) and _is_carboxyl_carbon(nb):
                    tail_atom = nb
                    break

        if head_atom is None and matched_atom_names:
            for atom in molecule.GetAtoms():
                if in_scope(atom) and atom.GetSymbol() == 'N':
                    head_atom = atom
                    break

        return head_atom, main_chain_atom, tail_atom

    head_atom, main_chain_atom, tail_atom = _resolve_mapping(use_match_filter)
    if matched_atom_names and (head_atom is None or main_chain_atom is None or tail_atom is None):
        head_fallback, main_fallback, tail_fallback = _resolve_mapping(False)
        if head_atom is None:
            head_atom = head_fallback
        if main_chain_atom is None:
            main_chain_atom = main_fallback
        if tail_atom is None:
            tail_atom = tail_fallback

    amber_mapping['HEAD_NAME'] = [_atom_name(head_atom)] if head_atom else []
    amber_mapping['TAIL_NAME'] = [_atom_name(tail_atom)] if tail_atom else []
    amber_mapping['MAIN_CHAIN'] = [_atom_name(main_chain_atom)] if main_chain_atom else []

    # 开始获取OMIT_NAME, 
    # 1. 与TAIL_NAME相邻的OH原子，包括O原子及其相邻的H原子
    for atom in molecule.GetAtoms():
        if _atom_name(atom) in amber_mapping['TAIL_NAME']:
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 1:
                    amber_mapping['OMIT_NAME'].append(_atom_name(n))
                    amber_mapping['OMIT_NAME'].extend([_atom_name(n) for n in n.GetNeighbors() if n.GetSymbol() == 'H'])

                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 0:
                    O_mapping = _atom_name(n) # O_mapping是O原子的名称
                    O_mappings.append(O_mapping)

    if amber_mapping['HEAD_NAME']:
        for atom in molecule.GetAtoms():
            if _atom_name(atom) == amber_mapping['HEAD_NAME'][0]:
                # 2. 与HEAD_NAME相邻的NH原子，有2个H，但是只删除一个H
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == 'H' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'N']) == 1:
                        amber_mapping['OMIT_NAME'].append(_atom_name(n))
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

def modify_prepin(prepin, atom_map, O_mappings, head_h_names=None):
    '''
    Modifies the prepin file using direct column updates.
    '''

    def _update_line(line, replacements):
        if len(line) < 18:
            return line
        if not line[:4].strip().isdigit():
            return line
        atom_name = line[6:10].strip()
        if atom_name not in replacements:
            return line
        new_name, new_type = replacements[atom_name]
        return f"{line[:6]}{new_name:<4}{line[10:12]}{new_type:<4}{line[16:]}"

    head_names = atom_map['HEAD_NAME']
    tail_names = atom_map['TAIL_NAME']
    main_chains = atom_map['MAIN_CHAIN']

    replacements = {}
    for head_name in head_names:
        replacements[head_name] = ('N', 'N')
    for tail_name in tail_names:
        replacements[tail_name] = ('C', 'C')
    for main_chain in main_chains:
        replacements[main_chain] = ('CA', 'CT')
    for O_mapping in O_mappings:
        replacements[O_mapping] = ('O', 'O')
    if head_h_names:
        for h_name in head_h_names:
            replacements[h_name] = (h_name, 'H')

    required = set(head_names) | set(tail_names) | set(main_chains)
    required.update(O_mappings)
    updated = {name: False for name in replacements}
    with open(prepin, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        atom_name = line[6:10].strip() if len(line) >= 10 else ''
        if atom_name in replacements:
            updated[atom_name] = True
        new_lines.append(_update_line(line, replacements))

    missing = [name for name in required if not updated.get(name, False)]
    if missing:
        raise ValueError(f'Failed to update atom names/types in prepin: {", ".join(missing)}')

    with open(prepin, 'w') as f:
        f.writelines(new_lines)

def write_amber_mc(amber_map, output_file):
    """
    Writes the Amber mapping to an mc file.
    """
    if not amber_map["HEAD_NAME"] or not amber_map["TAIL_NAME"] or not amber_map["MAIN_CHAIN"]:
        raise ValueError('HEAD_NAME/TAIL_NAME/MAIN_CHAIN not found; check SMILES mapping or pattern matching.')
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

def _format_net_charge(charge):
    if charge is None:
        return None
    try:
        charge_val = float(charge)
    except (TypeError, ValueError):
        raise ValueError(f'Invalid charge value: {charge}')
    if not charge_val.is_integer():
        raise ValueError(f'Net charge must be an integer for antechamber: {charge}')
    return str(int(charge_val))

def _ac_has_du(ac_file):
    try:
        with open(ac_file, 'r') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                parts = line.split()
                if parts and parts[-1] == 'DU':
                    return True
    except FileNotFoundError:
        return False
    return False

def _find_head_hydrogens(ac_file, head_names):
    if not head_names:
        return []
    name_to_index = {}
    index_to_name = {}
    index_to_type = {}
    bonds = {}
    try:
        with open(ac_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    idx = int(parts[1])
                    name = parts[2]
                    atype = parts[-1] if parts else ''
                    name_to_index[name] = idx
                    index_to_name[idx] = name
                    index_to_type[idx] = atype
                elif line.startswith('BOND'):
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    idx1 = int(parts[2])
                    idx2 = int(parts[3])
                    bonds.setdefault(idx1, set()).add(idx2)
                    bonds.setdefault(idx2, set()).add(idx1)
    except FileNotFoundError:
        return []

    h_names = set()
    for head_name in head_names:
        head_idx = name_to_index.get(head_name)
        if head_idx is None:
            continue
        for nb in bonds.get(head_idx, []):
            nb_name = index_to_name.get(nb)
            nb_type = index_to_type.get(nb, '')
            if nb_name and (nb_name.startswith('H') or (nb_type and nb_type[0].lower() == 'h')):
                h_names.add(nb_name)
    return sorted(h_names)

def sdf_to_ac(sdf, name=None, charge=None, atom_type='amber'):
    """将SDF转换为Amber mc文件"""
    '''
    antechamber -fi sdf -i R1A.sdf -bk R1A -fo ac -o R1A.ac -c bcc -at amber
    '''
    net_charge = _format_net_charge(charge)
    cmd = f'antechamber -fi sdf -i {sdf} -bk {name} -fo ac -o {name}.ac -c bcc -at {atom_type}'
    if net_charge is not None:
        cmd = f'{cmd} -nc {net_charge}'
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
    ac1 = sdf_to_ac(sdf1, args.name1, args.charge1, atom_type='amber')
    if _ac_has_du(ac1):
        print('Info: detected DU atom type with -at amber; retrying with -at gaff2.')
        ac1 = sdf_to_ac(sdf1, args.name1, args.charge1, atom_type='gaff2')
    ac2 = sdf_to_ac(sdf2, args.name2, args.charge2, atom_type='amber')
    if _ac_has_du(ac2):
        print('Info: detected DU atom type with -at amber; retrying with -at gaff2.')
        ac2 = sdf_to_ac(sdf2, args.name2, args.charge2, atom_type='gaff2')

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
    head_h_names_1 = _find_head_hydrogens(ac1, amber_map_1['HEAD_NAME'])
    head_h_names_2 = _find_head_hydrogens(ac2, amber_map_2['HEAD_NAME'])
    modify_prepin(prepin1, amber_map_1, O_mapping_1, head_h_names=head_h_names_1)
    modify_prepin(prepin2, amber_map_2, O_mapping_2, head_h_names=head_h_names_2)
    prep_to_frcmod(prepin1, f'frcmod.{args.name1}')
    prep_to_frcmod(prepin2, f'frcmod.{args.name2}')

    # combine the smi1 and smi2
    combined_charges = float(args.charge1) + float(args.charge2)
    combined_smiles = connect_molecules(args.smiles1, args.smiles2, f'{args.name1}_{args.name2}', nbonds=args.nbonds) # 生成连接后的SMILES
    combined_sdf = smi_to_sdf(process_smiles(combined_smiles), f'{args.name1}_{args.name2}') # 生成连接后的SDF
    combined_ac = sdf_to_ac(combined_sdf, f'{args.name1}_{args.name2}', combined_charges, atom_type='amber') # 生成连接后的AC
    if _ac_has_du(combined_ac):
        print('Info: detected DU atom type with -at amber; retrying with -at gaff2.')
        combined_ac = sdf_to_ac(combined_sdf, f'{args.name1}_{args.name2}', combined_charges, atom_type='gaff2') # 生成连接后的AC
    combined_pdb = ac_to_pdb(combined_ac, f'{args.name1}_{args.name2}') # 生成连接后的PDB
    combined_mol = read_molecule_from_pdb(combined_pdb)
    combined_matched_atoms = find_smiles_pattern(combined_pdb)
    combined_amber_map, combined_O_mapping = map_to_amber_mc(combined_mol, combined_matched_atoms, charge=combined_charges)
    modify_ac(combined_ac, combined_amber_map, combined_O_mapping)
    ac_to_frcmod(combined_ac, f'frcmod.{args.name1}_{args.name2}')

    # Change back to the original directory
    os.chdir(pwd)

    ## Example usage:
    # python generate_covalent_template.py --smiles1 "[*]C(C)C[C@@H](N)C(=O)O" --name1 R1A --smiles2 "[*]CC(C)C[C@@H](N)C(=O)O" --name2 R2A --output R1A_R2A --nbonds 2
