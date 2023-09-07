import os
import uuid
import shutil
import logging
from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor

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

    def generate_3d_structure_from_template(self, seq, output_pdb, template_pdb: str):
        '''
            Generate a 3D structure of a peptide from a template using Modeller.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                template_pdb (str): The path to the template PDB file.

            Returns:
                str: The path to the generated PDB file.
        '''
        # get absolute path of template pdb file
        template_pdb = os.path.abspath(template_pdb)
        pp = PrepareProt(seq, self.tmp_dir, method='modeller', template_pdb_file_path=template_pdb)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def de_novo_3d_structure(self, seq, output_pdb):
        '''
            Generate a de novo 3D structure of a peptide using ESMFold.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.

            Returns:
                str: The path to the generated PDB file.
        '''
        pp = PrepareProt(seq, self.tmp_dir, method='alphafold')
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def generate_3d_structure_from_sequence(self, seq, output_pdb):
        '''
            Generate a 3D structure of a peptide using Ambertools.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.

            Returns:
                str: The path to the generated PDB file.

            Note:
                This method is not recommended as the generated structure is not stable.
        '''
        pp = PrepareProt(seq, self.tmp_dir, method=None)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb
        
    
if __name__ == '__main__':
    seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'
    st = Structure(verbose=True)
    st.generate_3d_structure_from_template(seq=seq, 
                                           output_pdb='homology_model.pdb', 
                                           template_pdb='example/data/template.pdb')
    st.generate_3d_structure_from_sequence(seq=seq, 
                                           output_pdb='sequence.pdb')
    st.de_novo_3d_structure(seq=seq, 
                            output_pdb='denovo.pdb')