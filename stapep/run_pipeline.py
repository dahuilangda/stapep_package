import os
import argparse
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor, ProtParamsSeq

def md(args):
    """Run molecular dynamics simulation"""

    # Prepare the protein for molecular dynamics simulation
    pp = PrepareProt(args.seq, args.output, method=args.method, template_pdb_file_path=args.template_pdb)
    pp._gen_prmtop_and_inpcrd_file()

    # Run molecular dynamics simulation
    sim = Simulation(args.output)
    sim.setup(
            type=args.type, # 'explicit' or 'implicit'
            solvent=args.solvent,
            temperature=args.temperature, # K
            friction=args.friction, # ps^-1
            timestep=args.timestep, # ps
            interval=args.interval, # ps
            nsteps=args.nsteps
        )
    sim.minimize()
    sim.run()

def physicochemical(args):
    """Calculate physicochemical properties of the protein"""

    if args.type == 'explicit':
        topology_file = os.path.join(args.output, 'pep.prmtop')
    elif args.type == 'implicit':
        topology_file = os.path.join(args.output, 'pep_vac.prmtop')
    else:
        raise ValueError('Invalid type of simulation, must be either explicit or implicit')
    pcp = PhysicochemicalPredictor(sequence=args.seq, 
                                    topology_file=topology_file,
                                    trajectory_file=os.path.join(args.output, 'traj.dcd'),
                                    start_frame=args.start_frame)

    # save the mean structure of the trajectory
    mean_structure = 'mean_structure.pdb'
    pcp._save_mean_structure(os.path.join(args.output, mean_structure))

    return {
        'mean_bfactor': pcp.calc_mean_bfactor(),
        'mean_molsurf': pcp.calc_mean_molsurf(),
        'mean_gyrate': pcp.calc_mean_gyrate(),
        'hydrophobic_index': pcp.calc_hydrophobic_index(
            os.path.join(args.output, mean_structure)
        ),
        'psa': pcp.calc_psa(os.path.join(args.output, mean_structure)),
        'num_hbonds': pcp.calc_n_hbonds(),
        'helix_percent': pcp.calc_helix_percent(),
        'extend_percent': pcp.calc_extend_percent(),
        'loop_percent': pcp.calc_loop_percent(),
    }

def sequence(args):
    """Calculate sequence-based properties of the protein"""

    pps = ProtParamsSeq(seq=args.seq)

    return {
        'length': pps.seq_length,
        'weight': pps.weight,
        'hydrophobicity index': pps.hydrophobicity_index,
        'charge': pps.calc_charge(pH=args.ph),
        'charge_density': pps.calc_charge_density(pH=args.ph),
        'aromaticity': pps.aromaticity,
        'fraction_arginine': pps.fraction_arginine,
        'fraction_lysine': pps.fraction_lysine,
        'lyticity index': pps.lyticity_index,
        'isoelectric_point': pps.isoelectric_point
    }

def main(args):
    md(args)
    feats = physicochemical(args)
    feats.update(sequence(args))

    if args.permeability:
        import pickle
        feature_cols = ['length', 'weight', 'helix_percent',
            'loop_percent', 'mean_bfactor', 'mean_gyrate', 'hydrophobic_index',
            'num_hbonds', 'charge', 'aromaticity', 'isoelectric_point',
            'fraction_arginine', 'fraction_lysine', 'psa']
        
        X = pd.DataFrame(feats, index=[0])[feature_cols]
        model_path = f"{os.path.dirname(os.path.realpath(__file__))}/models/lgb_model.sav"
        model = pickle.load(open(model_path, 'rb'))
        feats['permeability'] = model.predict(X)[0]
        print(f"Predicted permeability: {feats['permeability']}")

    return feats

def get_args():
    parser = argparse.ArgumentParser(description='Run molecular dynamics simulation and extract features...')
    parser.add_argument('--seq', type=str, required=True, help='Sequence of the protein')
    parser.add_argument('--output', type=str, required=True, help='Output directory')
    parser.add_argument('--type', type=str, default='implicit', help='Type of simulation')
    parser.add_argument('--solvent', type=str, default='water', help='Solvent')
    parser.add_argument('--temperature', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--friction', type=float, default=1, help='Friction (ps^-1)')
    parser.add_argument('--timestep', type=float, default=0.002, help='Timestep (ps)')
    parser.add_argument('--interval', type=int, default=1000, help='Interval (ps)')
    parser.add_argument('--nsteps', type=int, default=5000000, help='Number of steps')
    parser.add_argument('--method', type=str, default='alphafold', help='Method to generate 3D structure, including alphafold, modeller or None')
    parser.add_argument('--template_pdb', type=str, default=None, help='Template PDB file for modeller, only used when method is modeller')
    parser.add_argument('--ph', type=float, default=7.0, help='pH')
    parser.add_argument('--start_frame', type=int, default=0, help='Start frame')
    parser.add_argument('--permeability', action='store_true', help='Calculate permeability using build-in model')

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    feats = main(args)
    df = pd.DataFrame(feats, index=[0])
    df.to_csv(os.path.join(args.output, 'feats.csv'), index=False)
    # print(df)
    print(f'Features saved to {os.path.join(args.output, "feats.csv")}')
