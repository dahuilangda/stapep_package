import os
import uuid
import shutil
import logging

import pandas as pd
from Bio import AlignIO
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist
import Bio.Align.substitution_matrices as matlist
try:
    from Bio.Data.PDBData import protein_letters_3to1 as aa3to1
except ImportError:
    print('Warning: Bio.Data.PDBData is deprecated, please use Bio.Data.SCOPData instead')
    from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor, SeqPreProcessing

class Structure(object):
    """
    The Structure class represents a structure and provides methods for generating 3D structures of peptides.

    Args:
        solvent (str, optional): The solvent to use for the simulation. Defaults to 'water'.
        save_tmp_dir (bool, optional): Whether to save the temporary directory used for the simulation. Defaults to False.
        verbose (bool, optional): Whether to enable verbose logging. Defaults to False.

    Examples:
        ```python
        # Create an instance of the Structure class
        structure = Structure()

        # Show the available solvent options
        structure.show_solvent_options()

        # Generate a 3D structure from a template
        structure.generate_3d_structure_from_template('Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ', 'output.pdb', 'template.pdb')

        # Generate a de novo 3D structure
        structure.de_novo_3d_structure('Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ', 'output.pdb')

        # Generate a 3D structure from a sequence
        structure.generate_3d_structure_from_sequence('ACDEFG', 'output.pdb')
        ```
    """
    def __init__(self, 
                 solvent: str='water',
                 save_tmp_dir: bool=False,
                 verbose: bool=False):
        self.solvent = solvent
        self.tmp_dir = os.path.join('/tmp', str(uuid.uuid4()))
        self.save_tmp_dir = save_tmp_dir
        self.verbose = verbose

        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        if verbose:
            logging.basicConfig(level=logging.INFO)
            if self.save_tmp_dir:
                logging.info(f'Temporarily directory: {self.tmp_dir}')

    def show_solvent_options(self):
        print('You can choose the following solvent options:')
        print('- water: default')
        print('- chloroform')
        print('- DMF: dimethylformamide')
        print('- DMSO: dimethyl sulfoxide')
        print('- ethanol')
        print('- acetone')

    def _check_solvent(self, solvent: str):
        if solvent not in ['water', 'chloroform', 'DMF', 'DMSO', 'ethanol', 'acetone']:
            raise ValueError(f'{solvent} is not a valid solvent option, please choose from the following options: water, chloroform, DMF, DMSO, ethanol, acetone')

    def _short_time_simulation(self, nsteps: int=100000):
        sim = Simulation(self.tmp_dir)
        sim.setup(type='implicit', 
                  solvent=self.solvent, 
                  temperature=300, 
                  friction=1, 
                  timestep=2, 
                  interval=10, 
                  nsteps=nsteps)
        if self.verbose:
            logging.info(f'Running short time simulation for {nsteps} steps')
        sim.minimize()
        sim.run()

    def _get_opt_structure(self, seq, pdb):
        pcp = PhysicochemicalPredictor(sequence=seq, 
                                        topology_file=os.path.join(self.tmp_dir, 'pep_vac.prmtop'), # topology file　(default: pep_vac.prmtop in the data folder)
                                        trajectory_file=os.path.join(self.tmp_dir, 'traj.dcd'), # trajectory file (default: traj.dcd in the data folder)
                                        start_frame=0) # start frame (default: 500)
        pcp._save_mean_structure(pdb)

    def _del_tmp_dir(self):
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

    def generate_3d_structure_from_template(self, 
                                            seq: str, 
                                            output_pdb: str, 
                                            template_pdb: str,
                                            additional_residues: dict=None):
        '''
            Generate a 3D structure of a peptide from a template using Modeller.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                template_pdb (str): The path to the template PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.
        '''
        spp = SeqPreProcessing(additional_residues=additional_residues)
        spp.check_seq_validation(seq)
        # get absolute path of template pdb file
        template_pdb = os.path.abspath(template_pdb)
        pp = PrepareProt(seq, self.tmp_dir, method='modeller', template_pdb_file_path=template_pdb)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def de_novo_3d_structure(self, 
                             seq: str, 
                             output_pdb: str,
                             additional_residues: dict=None,
                             proxy=None):
        '''
            Generate a de novo 3D structure of a peptide using ESMFold.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.
        '''
        spp = SeqPreProcessing(additional_residues=additional_residues)
        spp.check_seq_validation(seq)
        pp = PrepareProt(seq, self.tmp_dir, method='alphafold', additional_residues=additional_residues)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def generate_3d_structure_from_sequence(self, 
                                            seq: str, 
                                            output_pdb: str,
                                            additional_residues: dict=None):
        '''
            Generate a 3D structure of a peptide using Ambertools.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.

            Note:
                This method is not recommended as the generated structure is not stable.
        '''
        spp = SeqPreProcessing(additional_residues=additional_residues)
        spp.check_seq_validation(seq)
        pp = PrepareProt(seq, self.tmp_dir, method=None)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb
        
class AlignStructure(object):

    @staticmethod
    def convert_pdb_to_seq(id: str, pdb_file: str) -> list[str]:
        '''
            Convert a pdb file to a fasta file.
        '''
        res_list = AlignStructure._get_pdb_sequence(id, pdb_file)
        seq = [res[1] for res in res_list]
        seq = ''.join(seq)
        return seq

    @staticmethod
    def _get_pdb_sequence(id, pdb_file) -> list[tuple]:
        '''
            Return a list of tuples (idx, sequence).
            eg:[(6, 'P'),
                (7, 'D'),
                (8, 'I'),
                (9, 'F'),]
        '''
        parser = PDBParser()
        structure = parser.get_structure(id, pdb_file)
        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        return [_aainfo(r) for r in structure.get_residues() if is_aa(r)]

    @staticmethod
    def align(ref_pdb: str, pdb: str, output_pdb: str):
        '''
            Align structures using BioPython.

            Args:
                ref_pdb (str): The path to the reference PDB file.
                pdb (str): The path to the PDB file to align.
                output_pdb (str): The path to save the output PDB file.

            Returns:
                str: The path to the generated PDB file.
        '''

        try:
            import pymol
            from pymol import cmd
            cmd.reinitialize()
        except Exception as e:
            raise ImportError('Please install PyMOL to use this method. mamba install -c conda-forge pymol-open-source')

        # Load the PDB files
        pymol.cmd.load(ref_pdb, 'ref')
        pymol.cmd.load(pdb, 'denovo')

        # Perform alignment on alpha carbons (CA atoms)
        out = pymol.cmd.align('denovo and name CA', 'ref and name CA')
        rmsd, n_atoms, n_cycles, n_rmsd_pre, n_atom_pre, score, n_res = out

        # Make sure to update coordinates of the denovo structure
        pymol.cmd.alter_state(1, 'denovo', 'x, y, z = x, y, z')
        # Apply the transformation matrix after alignment
        pymol.cmd.matrix_copy('denovo', 'ref')
        # Save the aligned denovo structure to a new PDB file
        pymol.cmd.save(output_pdb, 'denovo')
        return rmsd

    # @staticmethod
    # def align(ref_pdb: str, pdb: str, output_pdb: str):
    #     '''
    #         Align structures using BioPython.

    #         Args:
    #             ref_pdb (str): The path to the reference PDB file.
    #             pdb (str): The path to the PDB file to align.
    #             output_pdb (str): The path to save the output PDB file.

    #         Returns:
    #             str: The path to the generated PDB file.
    #     '''

    #     ref_id = os.path.basename(ref_pdb).split('.')[0]
    #     pdb_id = os.path.basename(pdb).split('.')[0]


    #     parser = PDBParser()
    #     ref_structure = parser.get_structure(ref_id, ref_pdb)
    #     pdb_structure = parser.get_structure(pdb_id, pdb)
    #     ref_model = list(ref_structure.get_models())[0]
    #     pdb_model = list(pdb_structure.get_models())[0]

    #     resseq_A = AlignStructure._get_pdb_sequence(ref_id, ref_pdb)
    #     resseq_B = AlignStructure._get_pdb_sequence(pdb_id, pdb)
    #     sequence_A = AlignStructure.convert_pdb_to_seq(ref_id, ref_pdb)
    #     sequence_B = AlignStructure.convert_pdb_to_seq(pdb_id, pdb)

    #     alns = pairwise2.align.globalds(sequence_A, sequence_B, matlist.load("BLOSUM62"), -10.0, -0.5,
    #                                         penalize_end_gaps=(False, False) )
    #     best_aln = alns[0]
    #     aligned_A, aligned_B, score, begin, end = best_aln
    #     mapping = {}
    #     aa_i_A, aa_i_B = 0, 0
    #     for aa_aln_A, aa_aln_B in zip(aligned_A, aligned_B):
    #         if aa_aln_A == '-':
    #             if aa_aln_B != '-':
    #                 aa_i_B += 1
    #         elif aa_aln_B == '-':
    #             if aa_aln_A != '-':
    #                 aa_i_A += 1
    #         else:
    #             assert resseq_A[aa_i_A][1] == aa_aln_A
    #             assert resseq_B[aa_i_B][1] == aa_aln_B
    #             mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
    #             aa_i_A += 1
    #             aa_i_B += 1

    #     # Extract CA atoms using the helper function
    #     refe_ca_list = AlignStructure.get_CA_atoms_from_model(ref_model, list(mapping.keys()))
    #     mobi_ca_list = AlignStructure.get_CA_atoms_from_model(pdb_model, list(mapping.values()))

    #     # Superimpose matching residues
    #     try:
    #         si = Superimposer()
    #         si.set_atoms(refe_ca_list, mobi_ca_list)
    #         si.apply(pdb_structure.get_atoms())

    #         io = PDBIO()
    #         io.set_structure(pdb_structure)
    #         io.save(output_pdb)
    #         return output_pdb
    #     except Exception as e:
    #         print(e)

    @staticmethod
    def rmsd(ref_pdb: str, pdb: str):
        try:
            import pymol
            from pymol import cmd
            cmd.reinitialize()
        except Exception as e:
            print(e)
            return None

        pymol.cmd.load(ref_pdb, 'ref')
        pymol.cmd.load(pdb, 'denovo')
        out = pymol.cmd.align('ref and name CA', 'denovo and name CA')
        rmsd, n_atoms, n_cyles, n_rmsd_pre, n_atom_pre, score, n_res = out
        return rmsd

    @staticmethod
    def get_CA_atoms_from_model(model, residue_numbers):
        """
        Get CA atoms from a given model based on specified residue numbers.

        Args:
        - model: BioPython Model object
        - residue_numbers: List of residue numbers to extract CA atoms for

        Returns:
        - List of CA atoms
        """
        ca_atoms = []

        for chain in model:
            ca_atoms.extend(
                res['CA']
                for res in chain
                if res.id[1] in residue_numbers and 'CA' in res
            )
        return ca_atoms

if __name__ == '__main__':
    seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'
    st = Structure(verbose=True)
    st.generate_3d_structure_from_template(seq=seq, 
                                           output_pdb='example/data/homology_model.pdb', 
                                           template_pdb='example/data/template.pdb')
    
    AlignStructure.align(ref_pdb='example/data/template.pdb',
             pdb='example/data/homology_model.pdb',
             output_pdb='example/data/aligned.pdb')

    st.generate_3d_structure_from_sequence(seq=seq, 
                                           output_pdb='example/data/sequence.pdb')
    st.de_novo_3d_structure(seq=seq, 
                            output_pdb='example/data/denovo.pdb')