amino_acid_dict = {'CYS': 'C',  'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', # standard amino acids
                   'ILE': 'I',  'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', # standard amino acids
                   'GLY': 'G',  'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', # standard amino acids
                   'ALA': 'A',  'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', # standard amino acids
                   'PS3': 'S3', 'PS5': 'S5', 'PS8': 'S8', # stapled amino acids
                   'PR3': 'R3', 'PR5': 'R5', 'PR8': 'R8', # stapled amino acids
                   'ACE': 'Ac', 'NME': 'NH2',  # N/C terminal
                   'NLE': 'B', 'AIB': 'Aib', # non-standard amino acids
                   }

reversed_amino_acid_dict = {a: b for b, a in amino_acid_dict.items()}

amino_acid_smiles_dict = {'TRP': 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N',
                          'PHE': 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
                          'TYR': 'C1=CC(=CC=C1CC(C(=O)O)N)O',
                          'ALA': 'C[C@@H](C(=O)O)N',
                          'VAL': 'CC(C)[C@@H](C(=O)O)N',
                          'LEU': 'CC(C)C[C@@H](C(=O)O)N',
                          'ILE': 'CC[C@H](C)[C@@H](C(=O)O)N',
                          'MET': 'CSCC[C@@H](C(=O)O)N',
                          'PRO': 'C1C[C@H](NC1)C(=O)O',
                          'THR': 'C[C@H]([C@@H](C(=O)O)N)O',
                          'SER': 'C([C@@H](C(=O)O)N)O',
                          'CYS': 'C([C@@H](C(=O)O)N)S',
                          'GLN': 'C(CC(=O)N)[C@@H](C(=O)O)N',
                          'GLU': 'C(CC(=O)O)[C@@H](C(=O)O)N',
                          'ASN': 'C([C@@H](C(=O)O)N)C(=O)N',
                          'ASP': 'C([C@@H](C(=O)O)N)C(=O)O',
                          'LYS': 'C(CCN)C[C@@H](C(=O)O)N',
                          'ARG': 'C(C[C@@H](C(=O)O)N)CN=C(N)N',
                          'HIS': 'C1=C(NC=N1)C[C@@H](C(=O)O)N',
                          'GLY': 'C(C(=O)O)N',
                          'AIB': 'CC(C)(C(=O)O)N',
                          'NLE': 'CCCC[C@@H](C(=O)O)N',
                          'PS3': 'C=CC[C@](C)(N)C(=O)O',
                          'PS5': 'C=CCCC[C@](C)(N)C(=O)O',
                          'PS8': 'C=CCCCCCC[C@](C)(N)C(=O)O',
                          'PR3': 'C=CC[C@@](C)(N)C(=O)O',
                          'PR5': 'C=CCCC[C@@](C)(N)C(=O)O',
                          'PR8': 'C=CCCCCCC[C@@](C)(N)C(=O)O'}

li_dict = {'W': 32.606, 'F': 29.127, 
           'L': 23.398, 'I': 20.014, 
           'M': 15.854, 'Y': 14.647,
           'V': 13.285, 'P': 8.593,
           'C': 7.906,  'A': 3.207,
           'E': 0.392, 'T': 2.117,
           'D': -0.274, 'Q': 0.338,
           'S': -0.157,  'N': -0.523,
           'G': 0.0,  'R': 0.981,
           'H': 2.653,  'K': -3.774,
           'X': 34.32, 'S5': 34.32, 'S3': 12.05, 'S8': 65.555,
           'R3': 12.32, 'R5': 33.902, 'R8': 65.272,
           'Aib': 8.493, 'NLE': 24.442, 'B': 24.442}

hydrophilic_residues = ['Q','N','S','R','D','H','K','E']

weight_dict = {'W': 204.22899999999996,
               'F': 165.19199999999998,
               'Y': 181.19099999999997,
               'A': 89.09400000000001,
               'V': 117.148,
               'L': 131.17499999999995,
               'I': 131.17499999999995,
               'M': 149.215,
               'P': 115.13199999999999,
               'T': 119.11999999999999,
               'S': 105.09299999999999,
               'C': 121.16099999999999,
               'Q': 146.146,
               'E': 147.13,
               'N': 132.11899999999997,
               'D': 133.103,
               'K': 146.19,
               'R': 174.20399999999998,
               'H': 155.15699999999998,
               'G': 75.06700000000001,
               'Aib': 103.121,
               'B': 131.17499999999995,
               'S3': 129.159,
               'S5': 157.213,
               'X': 157.213,
               'S8': 199.29399999999998,
               'R3': 129.159,
               'R5': 157.213,
               'R8': 199.29399999999998} # with water
            
hydrophobic_dict = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
                    'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
                    'H': -3.2,'I': 4.5,  'L': 3.8,  'K': -3.9,
                    'M': 1.9, 'F': 2.8,  'P': -1.6, 'S': -0.8,
                    'T': -0.7,'W': -0.9, 'Y': -1.3, 'V': 4.2,
                    'X': 4.1, 'S5': 4.1, 'S3': 0.1, 'S8': 9.6, # staple residues
                    'R3': 0.2, 'R5': 4.0, 'R8': 9.6, # staple residues
                    'Aib': -0.5, 'B': 2.3} # non-standard residues