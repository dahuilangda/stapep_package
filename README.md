STAPEP: is a toolset for stapled peptides.
----------

Background
----------
Stapled peptides are peptide chains comprised of both natural and non-standard amino acids, possessing unique structures and functions. The process of stapling results in increased stability and may improve the permeability of select peptides.

Currently, there is no unified computational tool to handle stapled peptides. This project aims to provide a unified toolset, including: **building the 2D or 3D structure of stapled peptides**, **molecular dynamics**, and **extracting features** from the sequence and structure for other downstream tasks.

Installation
------------
    Prerequisites
        Install Anaconda
        Install CD-HIT: https://sites.google.com/view/cd-hit
        Clone this repository

    Install
        conda create -n stap python=3.9
        conda activate stap
        conda install -c conda-forge mamba
        mamba install -c conda-forge openmm parmed pandas numpy networkx biopython rdkit pymol-open-source scikit-learn lightgbm
        mamba install -c conda-forge -c ambermd ambertools

        python setup.py install

Usage
----------
### Sequence input rules
Supports all natural amino acids, 2 non-standard amino acids (Aib, NLE), and 6 stapled peptides (S3, S5, S8, R3, R5, R8).
* Use the single-letter amino acid abbreviation
* Single-letter amino acids can be written continuously or separated by '-'. For example: RRRRRR
* Non-single-letter amino acids must be separated by '-'. For example: Ac-BATPR8RRR-Aib-LBRR3FKRLQ-NH2
* The N-terminal can be started with 'Ac-', add Acetyl group
* The C-terminal can be ended with '-NH2', add Amide group
* Aib is 2-Aminoisobutyric acid
* B is the abbreviation of Nle, Nle is Norleucine
* R3, R5, R8 are right-handed stapled peptides
* S3, S5, S8 are left-handed stapled peptides
* X represents S5


### Generate structures from stapled peptide sequence
* To begin, 3D structures of natural peptides corresponding to the peptide sequences were generated utilizing [ESMFold](https://github.com/facebookresearch/esm)
* Next, the non-standard amino acids in the naturally occurring peptide structures need to be named correctly. This is a crucial step to ensure the accuracy and consistency of the peptide structures. The next step is to utilize tleap in AmberTools to finalize the structure of the stapled peptide.
* Ultimately, in order to obtain stable 3D structures of the stapled peptide, it is necessary to perform simulation using OpenMM. By conducting these simulations, a more comprehensive understanding of the peptide's properties and behaviors in a simulated environment can be obtained.


Perform molecular dynamics simulation to generate 3D structures of the stapled peptide, as shown in the following demo.
```python
from stapep.molecular_dynamics import PrepareProt, Simulation
seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
output = 'data' # Define the output directory

# Prepare the protein for molecular dynamics simulation
# method could be in ['alphafold', 'modeller', None]
# 'alphafold': use ESMFold to predict the structure
pp = PrepareProt(seq, output, method='alphafold')
# # 'modeller': use Modeller to predict the structure and use the template_pdb_file_path to specify the template file
# pp = PrepareProt(seq, output, method='modeller', template_pdb_file_path='demo_homology_model_template.pdb')
# # None: use the sequence to generate the structure directly by AmberTools
# pp = PrepareProt(seq, output, method=None)

pp._gen_prmtop_and_inpcrd_file() # Generate the prmtop and inpcrd files of Amber format
# Molecular dynamics simulation
sim = Simulation(output)
sim.setup(
        # 'explicit' or 'implicit'
        type='implicit', 
        # 'water' or 'chloroform'
        solvent='water', 
        # Kelvin
        temperature=300, 
        # ps^-1
        friction=1, 
        # fs
        timestep=2, 
        # ps (demo only)
        interval=100, 
        # 1 ns (demo only)
        nsteps=500000 
    )
sim.minimize()
sim.run()
```
_Demo: examples/ex1_molecular_dynamics.ipynb_


### Generate features from stapled peptide sequence
This tool enables the export of a variety of features based on peptide sequences, including length (excluding N-terminal and C-terminal), molecular weight, hydrophobicity index, charge, charge density, aromaticity, fraction arginine, fraction lysine, lyticity index, and isoelectric point. 

```python
from stapep.utils import ProtParamsSeq

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
pps = ProtParamsSeq(seq) # Initialize the class

# Get the features
print('length: ', pps.seq_length) # length of the peptide sequence (without N-terminal and C-terminal)
print('weight: ', pps.weight) # molecular weight
print('hydrophobicity index: ', pps.hydrophobicity_index) # hydrophobicity index (Kyte and Doolittle, 1982)
print('charge', pps.calc_charge(pH=7.0)) # charge at pH=7.0
print('charge_density', pps.calc_charge_density(pH=7.0)) # charge density at pH=7.0 (charge/length)
print('aromaticity', pps.aromaticity) # aromaticity (number of aromatic amino acids/length)
print('fraction_arginine', pps.fraction_arginine) # fraction of arginine (number of arginine/length)
print('fraction_lysine', pps.fraction_lysine) # fraction of lysine (number of lysine/length)
print('lyticity index: ', pps.lyticity_index) # lyticity index (Mourtada and Rida, 2019)
print('isoelectric_point: ', pps.isoelectric_point) # isoelectric point

# Output:
length:  19
weight:  2413.9530000000004
hydrophobicity index:  -0.2000000000000001
charge 5.996
charge_density 0.0024838926027143026
aromaticity 0.05263157894736842
fraction_arginine 0.2631578947368421
fraction_lysine 0.05263157894736842
lyticity index:  114.16199999999999
isoelectric_point:  11.999967765808105
```


The Lyticity Index is a valuable tool for discovering antimicrobial peptides. To broaden its utility, we expanded the index based on the groundwork laid by [Mouratada et al](https://doi.org/10.1038%2Fs41587-019-0222-z). (2019), to encompass six varieties of stapled peptides and two classes of non-standard amino acids in addition to S5. We subsequently re-measured the Lyticity Index values of all-natural peptides. The formula for calculating this index is as follows:
$$LI = \sum\limits_{i=1}^{n}(H_i + H_{i+4}) + \sum\limits_{i=1}^{n-1}(H_i + H_{i+3})$$
where $H_i$ is the hydrophobicity value of the i-th amino acid, and n is the number of amino acids.
| Amino Acid | $\Delta t_{\text{Gly}}$ | Amino Acid | $\Delta t_{\text{Gly}}$ |
| --- | --- | --- | --- |
| W | 32.606 | F | 29.127 |
| L | 23.398 | I | 20.014 |
| M | 15.854 | Y | 14.647 |
| V | 13.285 | P | 8.593 |
| C | 7.906 | A | 3.207 |
| E | 0.392 | T | 2.117 |
| D | -0.274 | Q | 0.338 |
| S | -0.157 | N | -0.523 |
| G | 0.0 | R | 0.981 |
| H | 2.653 | K | -3.774 |
| X (S5) | 34.32 | S3 | 12.05 |
| S8 | 65.555 | R3 | 12.32 |
| R5 | 33.902 | R8 | 65.272 |
| NLE | 24.442 | Aib | 8.493 |

$\Delta t_{\text{Gly}}$ is retention time relative to glycine at 25°C in reversed-phase chromatography (RP-HPLC).


In addition to calculating the lyticity index, we can also draw the image
```python
import os
pathname = 'data'
pps.plot_lyticity_index(os.path.join(pathname, 'lyticity_index.png')) # Draw the lyticity index image and save it to the specified path
```
![lyticity_index](https://github.com/dahuilangda/stapep_code/blob/master/example/img/Lyticity_index.png)


### Generate features from 3D structure of stapled peptide
The supported features include the average structure, as well as the Mean B-factor, molecular surface, Mean Gyration Radius, Hydrophobic Index, and 3D-PSA.
```python
from stapep.utils import PhysicochemicalPredictor

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
pathname = 'data' # The path where the topology file and trajectory file are located
pcp = PhysicochemicalPredictor(sequence=seq, 
                                topology_file=os.path.join(pathname, 'pep_vac.prmtop'), # topology file　(default: pep_vac.prmtop in the data folder)
                                trajectory_file=os.path.join(pathname, 'traj.dcd'), # trajectory file (default: traj.dcd in the data folder)
                                start_frame=0) # start frame (default: 500)

# Get the features
print('helix percent: ', pcp.calc_helix_percent()) # Calculate the percentage of helix in the secondary structure
print('sheet percent: ', pcp.calc_extend_percent()) # Calculate the percentage of sheet in the secondary structure
print('loop percent: ', pcp.calc_loop_percent()) # Calculate the percentage of loop in the secondary structure

# save the mean structure of the trajectory
pcp._save_mean_structure(os.path.join(pathname, 'mean_structure.pdb'))

# calculate the Mean B-factor, Molecular Surface, Mean Gyration Radius, Hydrophobic Index, and 3D-PSA
print('mean bfactor: ', pcp.calc_mean_bfactor())
print('mol surf: ', pcp.calc_mean_molsurf())
print('mean gyrate: ', pcp.calc_mean_gyrate())
print('hydrophobic index: ', pcp.calc_hydrophobic_index(os.path.join(pathname, 'mean_structure.pdb')))
print('psa: ', pcp.calc_psa(os.path.join(pathname, 'mean_structure.pdb')))
print('total number of hydrogen bonds: ', pcp.calc_n_hbonds())

# Output:
helix percent:  0.36699
sheet percent:  0.0
loop percent:  0.63301
mean bfactor:  49.358555117984295
mol surf:  1860.3414306676775
mean gyrate:  9.613029600296342
hydrophobic index:  -0.016666666666666757
psa:  933.903076171875
total number of hydrogen bonds:  6
```


We can obtain the SMILES of the stapled peptide by the average structure.
```python
# extract 2D structure of the peptide
smiles = pcp.extract_2d_structure(os.path.join(pathname, 'mean_structure.pdb'))
print(smiles)

# Output:
'[H]O[C@]([H])(C([H])([H])[H])[C@@]([H])(C(=O)N1C([H])([H])C([H])([H])C([H])([H])[C@@]1([H])C(=O)N([H])C1(C([H])([H])[H])C(O)N([H])[C@@]([H])(C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H])C(=O)N([H])[C@@]([H])(C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H])C(=O)N([H])[C@@]([H])(C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H])C(=O)N([H])C(C([H])([H])[H])(C([H])([H])[H])C(O)N([H])C([H])(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C(O)N([H])[C@@]([H])(C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])C(=O)N([H])C([H])(C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C(O)N([H])[C@@]([H])(C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H])C(=O)N([H])C(C(O)N([H])[C@]([H])(C(=O)N([H])[C@]([H])(C(=O)N([H])[C@]([H])(C(=O)N([H])[C@]([H])(C(=O)N([H])[C@]([H])(C(=O)O)C([H])([H])C([H])([H])C(=O)N([H])[H])C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])[H])C([H])([H])c2c([H])c([H])c([H])c([H])c2[H])(C([H])([H])[H])C([H])([H])C([H])C([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C1([H])[H])N([H])C(=O)[C@@]([H])(N([H])C(O)C([H])(N([H])C(O)C([H])([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C([H])([H])[H]'

```
_Demo: examples/ex2_extract_features.ipynb_


Finally, we can use the following command to automatically obtain the features of the stapled peptide and predict the **cell permeability** through the built-in machine learning model.
```python
python -m stapep.run_pipeline --seq "Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ" \
                       --output example/data \
                       --method alphafold \
                       --permeability 
    args:
        --seq: peptide sequence
        --output: output directory
        --type: type of simulation(explicit or implicit) # default: implicit
        --solvent: solvent(water or chloroform) # default: water
        --temperature: temperature(K) # default: 300K
        --friction: friction(ps^-1) # default: 1ps^-1
        --timestep: timestep(ps) # default: 0.002ps
        --interval: interval(ps) # default: 1000ps
        --nsteps: number of steps # default: 5000000 (100ns)
        --alphafold: use ESMFold to predict the structure
        --permeability: predict the permeability using built-in machine learning model
        --ph: pH # default: 7.0
        --start_frame: start frame # default: 0
```
_Features will be saved to "example/data/feats.csv"_

### Apply features to machine learning using custom dataset
A dataset comprised of 585 peptides was collected for the purpose of classifying the permeability of stapled peptides. The dataset includes structural information and experimental data for each stapled peptide, including natural peptides. Features were extracted from each peptide and summarized in "example/datasets/Stapled-peptide_permeability.csv".

Different Filters Applied to Peptide Dataset
1.  **Duplicate peptide removal:** Duplicate peptides can arise from experimental errors, technical artifacts, or other factors, and their presence in the dataset can skew the results of downstream analyses. Removing duplicate peptides ensures that the dataset is not biased towards certain peptides and that the analyses are based on a representative set of unique peptides.
    
2.  **Length restriction:** Restricting the length of the peptides to a specific range can help to remove irrelevant peptides and focus on peptides that are more relevant to the research question. For example, if the research question is about cell penetration, then restricting the peptides to a certain length range can help to eliminate peptides that are too short or too long to penetrate cells.
    
3.  **Clustering-based filtering:** Clustering the peptides based on similarity can help to remove redundant peptides and identify representative peptides. This can be useful in cases where the dataset contains a large number of peptides with similar sequences or structures, as it can reduce the complexity of the dataset and improve the interpretability of the results.
```python
# Instantiate a Filter object by passing the input DataFrame and specifying the sequence and name columns.
import pandas as pd
from stapep.filter import Filter

df = pd.read_csv('stapep/example/datasets/Stapled-peptide_permeability.csv')

filter_ = Filter(df, seq_col='SEQUENCE', name_col='name')
print('Original', df.shape)

# Remove duplicate peptides.
df = filter_.drop_duplicate()
print('After drop duplicate', df.shape)

# Filter out peptide sequences based on their length (default: 10-50).
df = filter_.filter_seq_length()
print('After filter seq length', df.shape)

# # Filter class also provides a clustering method to filter peptide sequences by similarity threshold (CD-HIT program is required)
# df = filter_.filter_identity_threshold(0.9)
# print('After filter identity threshold', df.shape)

# Save the filtered dataset
df.to_csv('../datasets/Stapled-peptide_permeability_filtered.csv', index=False)

# Output:
Original (585, 21)
After drop duplicate (564, 21)
After filter seq length (506, 22)
```


Next, we use the following features to train a machine learning model to predict the permeability of stapled peptides
```python
# We selected these features for machine learning
feature_cols = ['length', 'weight', 'helix_percent',
       'loop_percent', 'mean_bfactor', 'mean_gyrate', 'hydrophobic_index',
       'num_hbonds', 'charge', 'aromaticity', 'isoelectric_point',
       'fraction_arginine', 'fraction_lysine', 'psa']

# drop rows with nan
df[feature_cols].dropna(inplace=True)

value_col = 'Permeability'
df['Value'] = df[value_col].replace({'Permeable': 1, 'Impermeable': 0})

# convert cols to number
for col in feature_cols:
       df[col] = pd.to_numeric(df[col], errors='coerce')

X_data = df[feature_cols]
y_data = np.array(list(df['Value'].values))

# np.array to dataframe
X = pd.DataFrame(X_data, columns=feature_cols)
y = pd.DataFrame(y_data, columns=['Value'])

X.isnull().sum().sort_values(ascending=False), y.isnull().sum().sort_values(ascending=False)

# split data into train and test
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 42)
```

We use GridSearchCV to find the optimal parameters, and then use the optimal parameters to train the model
```python
import lightgbm as lgb
from sklearn.model_selection import GridSearchCV

lgb_cls = lgb.LGBMClassifier(objective='binary', random_state=42)

param_grid = {
    'num_leaves': [15, 31, 50], # the maximum number of leaf nodes for each decision tree
    'learning_rate': [0.01, 0.05, 0.1], # the learning rate for each decision tree
    'n_estimators': [500, 1000, 2000], # how many decision trees to use
    'subsample': [0.6, 0.8, 1.0], # the proportion of samples randomly sampled during each decision tree training
    'colsample_bytree': [0.6, 0.8, 1.0] # the proportion of features randomly sampled during each decision tree training
}

grid_search = GridSearchCV(estimator=lgb_cls, param_grid=param_grid, cv=5)
grid_search.fit(X_train, y_train)

print("Best hyperparameters: ", grid_search.best_params_)
print("Training accuracy: ", grid_search.best_score_)

# Output:
Best hyperparameters:  {'colsample_bytree': 0.6, 'learning_rate': 0.1, 'n_estimators': 1000, 'num_leaves': 15, 'subsample': 0.6}
Training accuracy:  0.819074074074074
```

Here are the results of the model validation
```python
test_pred = grid_search.predict(X_test)
# ROC AUC
from sklearn.metrics import roc_auc_score
print('ROC AUC: ', roc_auc_score(y_test, test_pred))
# Accuracy
from sklearn.metrics import accuracy_score
print('Accuracy: ', accuracy_score(y_test, test_pred))
# Recall
from sklearn.metrics import recall_score
print('Recall: ', recall_score(y_test, test_pred))
# F1
from sklearn.metrics import f1_score
print('F1: ', f1_score(y_test, test_pred))

# Output:
ROC AUC:  0.8512572533849129
Accuracy:  0.8529411764705882
Recall:  0.8727272727272727
F1:  0.8648648648648648
```

Finally, we output feature importance and discuss the meaning of each feature
```python
import matplotlib.pyplot as plt
import seaborn as sns

# Get the feature importances from the best estimator
importance = grid_search.best_estimator_.feature_importances_

# Sort the indices of the feature importances in descending order
sorted_idx = importance.argsort()[::-1]

# Create a figure with a specified size
plt.figure(figsize=(10, 6))

# Create a bar plot of feature importances with a consistent color
sns.barplot(x=importance[sorted_idx], y=X_train.columns[sorted_idx], color='darkorange')

# Set plot title, labels, and adjust layout
plt.title('Feature Importance')
plt.xlabel('Importance Score')
plt.ylabel('Feature')
plt.tight_layout()

# Display the plot
plt.show()
```
![Feature Importance](https://github.com/dahuilangda/stapep_code/blob/master/example/img/Feature_importance.png)

### Factors Influencing Peptide Membrane Permeability
The feature importances suggest that the following factors are important for predicting the membrane permeability of peptides:

1.  **Hydrophobic Index:** This is a measure of the hydrophobicity or water-repelling nature of a peptide. Peptides with a higher hydrophobic index tend to have better membrane permeability, as they can more easily cross the hydrophobic barrier of the cell membrane.
    
2.  **Weight:** The molecular weight of a peptide can also influence its membrane permeability. Smaller peptides are generally more permeable than larger ones.
    
3.  **PSA (Polar Surface Area):** This is a measure of the surface area of a peptide that is exposed to polar solvents. Peptides with a lower PSA tend to have better membrane permeability, as they can more easily cross the non-polar, hydrophobic region of the cell membrane.
    
4.  **Mean Gyrate and Mean Bfactor:** These are measures of the flexibility and mobility of a peptide, respectively. Peptides that are more flexible and mobile may be better able to adopt the conformation required for membrane penetration.
    
5.  **Helix Percent and Loop Percent:** These are measures of the secondary structure of a peptide. Peptides that have a higher proportion of alpha helix structure and lower proportion of loop structure tend to have better membrane permeability.
    
6.  **Isoelectric Point:** This is the pH at which a peptide has a net charge of zero. Peptides that have a neutral charge at physiological pH (7.4) tend to have better membrane permeability.
    
7.  **Length, Aromaticity, Charge, Fraction Lysine, Fraction Arginine, and Num Hbonds:** These are other physicochemical properties of the peptide that may influence its membrane permeability. For example, peptides with a higher positive charge or higher fraction of basic amino acids such as lysine and arginine may have better membrane permeability due to interactions with the negatively charged cell membrane. Peptides with a higher number of hydrogen bonds may also have better membrane permeability due to stronger interactions with the cell membrane.