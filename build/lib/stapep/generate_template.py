import os 
import re
from rdkit import Chem
from rdkit.Chem import AllChem

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

def find_smiles_pattern(pdb, smiles="NCC(O)O"):
    """在PDB文件内容中寻找SMILES模式"""
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

    if matches:
        # 输出所有匹配的原子名称, 而不是原子索引，也不是原子类型
        pattern_with_map = Chem.MolFromSmiles(smiles)
        mapped_atom_indices = []
        if pattern_with_map:
            mapped_atom_indices = [
                idx for idx, atom in enumerate(pattern_with_map.GetAtoms())
                if atom.GetAtomMapNum() > 0
            ]
        use_atom_map = len(mapped_atom_indices) > 0
        atom_names = []
        for match in matches:
            if use_atom_map:
                atoms = [mol.GetAtomWithIdx(match[idx]) for idx in mapped_atom_indices]
                atom_names.extend([atom.GetPDBResidueInfo().GetName().strip() for atom in atoms])
                break
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            atom_names.extend([atom.GetPDBResidueInfo().GetName().strip() for atom in atoms])
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
        'HEAD_NAME': head_name,
        'TAIL_NAME': tail_name,
        'MAIN_CHAIN': main_chain,
        'OMIT_NAME': [] if omit_name is None else omit_name,
        'PRE_HEAD_TYPE': pre_head_type,
        'PRE_TAIL_TYPE': pre_tail_type,
        'CHARGE': charge
    }

    matched_atom_names = {name.strip() for name in matched_atoms} if matched_atoms else set()

    O_mapping = None

    use_match_filter = bool(matched_atom_names)
    for atom in molecule.GetAtoms():
        atom_name = atom.GetPDBResidueInfo().GetName().strip()
        if use_match_filter and atom_name not in matched_atom_names:
            continue

        # N原子作为头部：只连接一个重原子(通常是C)，且带有氢(显式或隐式)
        if atom.GetSymbol() == 'N':
            heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
            c_neighbors = [n for n in heavy_neighbors if n.GetSymbol() == 'C']
            if len(heavy_neighbors) == 1 and len(c_neighbors) == 1:
                amber_mapping['HEAD_NAME'] = atom_name

        # C原子: 相邻原子为=O，-O，C
        elif atom.GetSymbol() == 'C' and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) == 2 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']) == 1:
            amber_mapping['TAIL_NAME'] = atom_name

        # 链接C和N的原子为主链原子
        if atom.GetSymbol() == 'C' and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'N']) == 1 and \
            len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) == 0:
            amber_mapping['MAIN_CHAIN'] = atom_name

    if amber_mapping['HEAD_NAME'] is None and matched_atom_names:
        for atom in molecule.GetAtoms():
            atom_name = atom.GetPDBResidueInfo().GetName().strip()
            if atom_name in matched_atom_names and atom.GetSymbol() == 'N':
                amber_mapping['HEAD_NAME'] = atom_name
                break

    # 开始获取OMIT_NAME, 
    # 1. 与TAIL_NAME相邻的OH原子，包括O原子及其相邻的H原子
    for atom in molecule.GetAtoms():
        if atom.GetPDBResidueInfo().GetName().strip() == amber_mapping['TAIL_NAME']:
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 1:
                    amber_mapping['OMIT_NAME'].append(n.GetPDBResidueInfo().GetName().strip())
                    amber_mapping['OMIT_NAME'].extend([n.GetPDBResidueInfo().GetName().strip() for n in n.GetNeighbors() if n.GetSymbol() == 'H'])

                if n.GetSymbol() == 'O' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'H']) == 0:
                    O_mapping = n.GetPDBResidueInfo().GetName().strip()

    for atom in molecule.GetAtoms():
        if atom.GetPDBResidueInfo().GetName().strip() == amber_mapping['HEAD_NAME']:
            # 2. 与HEAD_NAME相邻的NH原子，有2个H，但是只删除一个H
            for n in atom.GetNeighbors():
                if n.GetSymbol() == 'H' and len([n for n in n.GetNeighbors() if n.GetSymbol() == 'N']) == 1:
                    amber_mapping['OMIT_NAME'].append(n.GetPDBResidueInfo().GetName().strip())
                    break

    return amber_mapping, O_mapping

def modify_prepin(prepin, atom_map, O_mapping):
    '''
    amber_map: {'HEAD_NAME': 'N3', 'TAIL_NAME': 'C6', 'MAIN_CHAIN': 'C5', 'OMIT_NAME': ['O2', 'H12', 'H10'], 'PRE_HEAD_TYPE': 'C', 'PRE_TAIL_TYPE': 'N', 'CHARGE': '0.0'}
    '''
    
    # 使用sed命令修改prepin文件
    # # 1. 修改N3为N
    head_name = atom_map['HEAD_NAME']
    if head_name is None:
        raise ValueError('HEAD_NAME is None; check SMILES mapping or pattern matching.')
    if len(head_name) == 2:
        head_name = f'{head_name} '
    cmd = f"sed -i 's/{head_name}   NT/N     N /g' {prepin}"
    print(cmd)
    os.system(cmd)

    # # 2. 修改C6为C
    tail_name = atom_map['TAIL_NAME']
    if tail_name is None:
        raise ValueError('TAIL_NAME is None; check SMILES mapping or pattern matching.')
    if len(tail_name) == 2:
        tail_name = f'{tail_name} '
    cmd = f"sed -i 's/{tail_name}   C/C     C/g' {prepin}"
    print(cmd)
    os.system(cmd)

    # # 3. 修改O1为O
    if O_mapping:
        if len(O_mapping) == 2:
            O_mapping = f'{O_mapping} '
        cmd = f"sed -i 's/{O_mapping}   O/O     O/g' {prepin}"
        print(cmd)
        os.system(cmd)

    # # 4. 修改C5位CA
    main_chain = atom_map['MAIN_CHAIN']
    if main_chain is None:
        raise ValueError('MAIN_CHAIN is None; check SMILES mapping or pattern matching.')
    if len(main_chain) == 2:
        main_chain = f'{main_chain} '
    cmd = f"sed -i 's/{main_chain}   CT/CA    CT/g' {prepin}"
    print(cmd)
    os.system(cmd)

def write_amber_mc(amber_map, output_file):
    """
    将Amber映射写入到mc文件中。
    Args:
    amber_map (dict): 包含HEAD_NAME, TAIL_NAME和MAIN_CHAIN字段的字典。
    output_file (str): 输出文件的路径。
    """
    with open(f'{output_file}.mc', 'w') as f:
        f.write(f'HEAD_NAME {amber_map["HEAD_NAME"]}\n')
        f.write(f'TAIL_NAME {amber_map["TAIL_NAME"]}\n')
        f.write(f'MAIN_CHAIN {amber_map["MAIN_CHAIN"]}\n')
        # f.write(f'OMIT_NAME {" ".join(amber_map["OMIT_NAME"])}\n')
        for omit_name in amber_map["OMIT_NAME"]:
            f.write(f'OMIT_NAME {omit_name}\n')
        f.write(f'PRE_HEAD_TYPE {amber_map["PRE_HEAD_TYPE"]}\n')
        f.write(f'PRE_TAIL_TYPE {amber_map["PRE_TAIL_TYPE"]}\n')
        f.write(f'CHARGE {amber_map["CHARGE"]}\n')
    return f'{output_file}.mc'

def write_std_amber_mc(amber_map, output_file):
    """
    将Amber映射写入到mc文件中。
    Args:
    amber_map (dict): 包含HEAD_NAME, TAIL_NAME和MAIN_CHAIN字段的字典。
    output_file (str): 输出文件的路径。
    """
    with open(f'{output_file}.mc', 'w') as f:
        f.write(f'HEAD_NAME N\n')
        f.write(f'TAIL_NAME C\n')
        f.write(f'MAIN_CHAIN CA\n')
        # f.write(f'OMIT_NAME {" ".join(amber_map["OMIT_NAME"])}\n')
        for omit_name in amber_map["OMIT_NAME"]:
            f.write(f'OMIT_NAME {omit_name}\n')
        f.write(f'PRE_HEAD_TYPE {amber_map["PRE_HEAD_TYPE"]}\n')
        f.write(f'PRE_TAIL_TYPE {amber_map["PRE_TAIL_TYPE"]}\n')
        f.write(f'CHARGE {amber_map["CHARGE"]}\n')
    return f'{output_file}.mc'

def smi_to_sdf(smi, name=None):
    """将SMILES转换为3D坐标"""
    '''
    obabel -ismi R1A.smi -osdf -O R1A.sdf --gen3D
    '''
    with open(f'{name}.smi', 'w') as f:
        f.write(_strip_atom_map(smi))

    cmd = f'obabel -ismi {name}.smi -osdf -O {name}.sdf --gen3D'
    print(cmd)
    os.system(cmd)
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

def sdf_to_ac(sdf, name=None, charge=None):
    """将SDF转换为Amber mc文件"""
    '''
    antechamber -fi sdf -i R1A.sdf -bk R1A -fo ac -o R1A.ac -c bcc -at amber
    '''
    net_charge = _format_net_charge(charge)
    cmd = f'antechamber -fi sdf -i {sdf} -bk {name} -fo ac -o {name}.ac -c bcc -at amber'
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
    """将Amber mc文件转换为prepin文件"""
    '''
    prepgen -i R1A.ac -o R1A.prepin -m R1A.mc -rn R1A
    '''
    if not name:
        name = mc.replace('.mc', '')
    cmd = f'prepgen -i {ac} -o {prepin} -m {mc} -rn {name}'
    print(cmd)
    os.system(cmd)
    return prepin

def prep_to_frcmod(prep, frcmod=None):
    """将prepin文件转换为frcmod文件"""
    '''
    parmchk2 -i R1A.prepin -f prepi -o frcmod.R1A -a Y
    '''
    cmd = f'parmchk2 -i {prep} -f prepi -o {frcmod} -a Y'
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate non-standard amino acid to Amber mc file')
    parser.add_argument('--smiles', help='PDB file')
    parser.add_argument('--name', help='Name of the non-standard amino acid, 3 characters long, For example: R1A')
    parser.add_argument('--charge', help='Charge of the non-standard amino acid', default='0.0')
    parser.add_argument('--output', help='Output directory', default=None)
    args = parser.parse_args()

    print(f'args.smiles: {args.smiles}')

    # 判断name是否为3个字符
    if len(args.name) != 3:
        raise ValueError('Name of the non-standard amino acid must be 3 characters')
    
    args.name = args.name.upper()

    if not args.output:
        args.output = args.name

    # 新建name文件夹，如果已存在则退出
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        raise ValueError(f'{args.name} already exists')
    
    # cd到name文件夹
    pwd = os.getcwd()
    os.chdir(args.output)

    sdf = smi_to_sdf(args.smiles, args.name)
    ac = sdf_to_ac(sdf, args.name, args.charge)
    pdb = ac_to_pdb(ac, args.name)

    pattern_smiles = args.smiles if _smiles_has_atom_map(args.smiles) else "NCC(O)O"
    matched_atoms = find_smiles_pattern(pdb, smiles=pattern_smiles)
    mol = read_molecule_from_pdb(pdb)
    amber_map, O_mapping = map_to_amber_mc(mol, matched_atoms, charge=args.charge)
    print(f'amber_map: {amber_map}', f'O_mapping: {O_mapping}')
    mc = write_amber_mc(amber_map, args.name)
    prepin = mc_to_prepin(ac, mc, f'{args.name}.prepin')
    modify_prepin(prepin, amber_map, O_mapping)
    prep_to_frcmod(prepin, f'frcmod.{args.name}')

    # cd到原来的目录
    os.chdir(pwd)
    
    # python generate_template.py --smiles "N=C(NCCC[C@H](N)C(=O)O)N[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O" --name R1A
