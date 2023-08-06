import os
import sys
sys.path.append('..')

import warnings
warnings.filterwarnings('ignore')

import pickle
import pandas as pd
from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor, ProtParamsSeq

import argparse


def main(seq, output):
    # Prepare the protein for molecular dynamics simulation
    if not os.path.exists(os.path.join(output, 'pep_vac.prmtop')):
        pp = PrepareProt(seq, output, alphafold=True)
        pp._gen_prmtop_and_inpcrd_file()

        sim = Simulation(output)
        sim.setup(
                type='implicit', # 'explicit' or 'implicit'
                solvent='water', # 'water' or 'chloroform'
                temperature=300, # Kelvin
                friction=1, # ps^-1
                timestep=2, # fs
                interval=1000,
                nsteps=5000000 # 10 ns
            )
        sim.minimize()
        sim.run()


    # Initialize the ProtParamsSeq class
    pps = ProtParamsSeq(seq)

    # Initialize the PhysicochemicalPredictor class, which automatically loads trajectory using pytraj.
    pcp = PhysicochemicalPredictor(sequence=seq, 
                                    topology_file=os.path.join(output, 'pep_vac.prmtop'),
                                    trajectory_file=os.path.join(output, 'traj.dcd'),
                                    start_frame=500)

    # Save the mean structure of the trajectory
    pcp._save_mean_structure(os.path.join(output, 'mean_structure.pdb'))

    # We selected these features for machine learning
    feature_cols = ['length', 'weight', 'helix_percent',
        'loop_percent', 'mean_bfactor', 'mean_gyrate', 'hydrophobic_index',
        'num_hbonds', 'charge', 'aromaticity', 'isoelectric_point',
        'fraction_arginine', 'fraction_lysine', 'psa']

    # calculate the features
    length = pps.seq_length
    weight = pps.weight
    helix_percent = pcp.calc_helix_percent()
    loop_percent = pcp.calc_loop_percent()
    mean_bfactor = pcp.calc_mean_bfactor()
    mean_gyrate = pcp.calc_mean_gyrate()
    hydrophobic_index = pcp.calc_hydrophobic_index(os.path.join(output, 'mean_structure.pdb'))
    num_hbonds = pcp.calc_n_hbonds()
    charge = pps.calc_charge(pH=7.0)
    aromaticity = pps.aromaticity
    isoelectric_point = pps.isoelectric_point
    fraction_arginine = pps.fraction_arginine
    fraction_lysine = pps.fraction_lysine
    psa = pcp.calc_psa(os.path.join(output, 'mean_structure.pdb'))

    # put the features into a DataFrame
    df = pd.DataFrame([[length, weight, helix_percent, loop_percent, 
                        mean_bfactor, mean_gyrate, hydrophobic_index, num_hbonds, 
                        charge, aromaticity, isoelectric_point, fraction_arginine, 
                        fraction_lysine, psa]], columns=feature_cols)
    

    # load the machine learning model
    model = pickle.load(open('../models/lgb_model.sav', 'rb'))

    # predict the permeability
    permeability_prob = model.predict_proba(df[feature_cols])[0][1]
    permeability_ = 'permeable' if permeability_prob > 0.5 else 'impermeable'
    print(f'The permeability of "{seq}" is "{permeability_} ({round(permeability_prob,3)})".')

    df['permeability'] = f'{permeability_} ({round(permeability_prob,3)})'
    df.to_csv(os.path.join(output, 'permeability.csv'), index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--seq', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)

    args = parser.parse_args()
    main(args.seq, args.output)