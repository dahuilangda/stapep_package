STAPEP: is a toolset for stapled peptides.
----------

Background
----------
Stapled peptides are peptide chains comprised of both natural and non-standard amino acids, possessing unique structures and functions. The process of stapling results in increased stability and may improve the permeability of select peptides.

Currently, there is no unified computational tool to handle stapled peptides. This project aims to provide a unified toolset, including: **building the 2D or 3D structure of stapled peptides**, **molecular dynamics**, and **extracting features** from the sequence and structure for other downstream tasks.

### News
- [2024/11/6] Our research article has been published in the Journal of Chemical Information and Modeling. For detailed insights, you can access the paper [here](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01718).

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
        pip install transformers
        pip install torch torchvision torchaudio

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
Our software offers a streamlined workflow for modeling and optimizing peptide structures, even those containing non-standard amino acids. Here's a step-by-step breakdown of the process:

1. **Initial Conversion**: Non-standard amino acids are initially converted to ALA (alanine) to facilitate the modeling of natural peptides.

2. **Peptide Structure Modeling**:
    - **Modeller - Homology Modeling**: Utilizes the principle that proteins with similar sequences often have analogous 3D structures. Using a known structure as a template, Modeller predicts the structure of the target protein.
    - **ESMFold - De Novo Prediction**: Without relying on a template, ESMFold predicts protein structures solely based on the peptide sequence.
    - **AmberTools - From Sequence to Structure**: Post peptide sequence definition, AmberTools' tleap module generates an initial 3D structure.

3. **Reconversion of Amino Acids**: Once the natural peptide structure is modeled, the names of the non-standard amino acids, previously converted to ALA, are reverted back to their original non-standard naming.

4. **Structural Completion with tleap**: Utilizing tleap in AmberTools, the peptide structure is finalized, ensuring the inclusion and correct positioning of non-standard amino acids.

5. **Dynamics Optimization with OpenMM**: To achieve stable and refined 3D structures of peptides, simulations using OpenMM are performed. This step provides a comprehensive understanding of the peptide's properties and behaviors within a simulated environment.

With this toolkit, users are equipped with the means to effortlessly model and optimize peptide structures, catering to both standard and non-standard amino acid scenarios, and setting the stage for advanced research applications.


Perform short time (100ps) molecular dynamics simulation to generate 3D structures of the stapled peptide, as shown in the following demo.

1. homology modeling using Modeller
```python
from stapep.structure import Structure, AlignStructure

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
st = Structure(verbose=True) # Initialize the structure class

# Generate the 3D structure of the stapled peptide using Modeller
st.generate_3d_structure_from_template(
        seq=seq,
        output_pdb='data/homology_model.pdb', 
        template_pdb='data/template.pdb')

# align structure to template (Optional)
AlignStructure.align(
        ref_pdb='data/template.pdb',
        pdb='data/homology_model.pdb',
        output_pdb='data/aligned.pdb')
```
2. de novo prediction using ESMFold
```python
from stapep.structure import Structure

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
st = Structure(verbose=True) # Initialize the structure class
st.de_novo_3d_structure(
        seq=seq, 
        output_pdb='data/de_novo.pdb')
```
3. Sequence to structure using AmberTools (Not recommended)
```python
from stapep.structure import Structure

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
st = Structure(verbose=True) # Initialize the structure class
st.generate_3d_structure_from_sequence(
        seq=seq, 
        output_pdb='data/from_sequence.pdb')
```
_Demo: examples/ex0_structure.ipynb_

#### Advanced Molecular Dynamics Simulation
While the structure module provides a quick approach to obtaining peptide structures through short-time molecular dynamics, users who wish for more accurate and stable features can opt for extended molecular dynamics simulations as described below:
```python
from stapep.molecular_dynamics import PrepareProt, Simulation
seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ' # Define the peptide sequence
output = 'data' # Define the output directory

# Prepare the protein for molecular dynamics simulation
# method could be in ['alphafold', 'modeller', None]
# Note: if you want to use the modeller method, you need to specify the template_pdb_file_path
pp = PrepareProt(seq, output, method='alphafold')
# pp = PrepareProt(seq, output, method='modeller', template_pdb_file_path='demo_homology_model_template.pdb')
# pp = PrepareProt(seq, output, method=None)

pp._gen_prmtop_and_inpcrd_file() # Generate the prmtop and inpcrd files of Amber format
# Molecular dynamics simulation
sim = Simulation(output)
sim.setup(type='implicit', # 'explicit' or 'implicit'
          solvent='water', # 'water' or 'chloroform'
          temperature=300, # Kelvin
          friction=1, # ps^-1
          timestep=2, # fs
          interval=100, # ps (demo only)
          nsteps=500000 # 1 ns (demo only)
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
lyticity index:  564.745
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
pps.plot_lyticity_index(os.path.join(pathname, 'lyticity_index.svg')) # Draw the lyticity index image and save it to the specified path
```
![lyticity_index](https://github.com/dahuilangda/stapep_package/blob/master/stapep/example/img/Lyticity_index.png)


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
        --method: method of structure modeling (alphafold, modeller, or None) # default: alphafold
        --template_pdb_file_path: template pdb file path # default: None, only used when method is modeller
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

### Experimental Feature: Handling Peptides with Novel Non-Standard Amino Acids
Taking glycosylated arginine peptides as an example, we provide a quick and simple method to obtain the parameters of non-standard amino acids. Based on the obtained parameters, subsequent operations such as 3D structure modeling, feature extraction, and machine learning can be performed.(Note: You should be careful to observe whether these parameters are reasonable in theory rather than using them without consideration.)

```bash
# Generate parameter of the glycosylated arginine peptide using the following command
python -m stapep.generate_template \
    --smiles "N=C(NCCC[C@H](N)C(=O)O)N[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O" \
    --name "R1A" \
    --output "./data/R1A" \
    --charge 0.0

args:
    --smiles: SMILES of the non-standard amino acid
    --name: name of the non-standard amino acid (Only 3 characters are allowed)
    --output: output directory
    --charge: charge of the non-standard amino acid (default: 0.0)
```

The SMILES structure of the non-standard amino acid should include the complete structure, including COOH and NH2 on the main chain, as well as all atoms on the side chain. As shown in the figure below.

![R1A](https://github.com/dahuilangda/stapep_package/blob/master/stapep/example/img/non_standard_aa.png)

After running, you can find the generated parameter files in the specified directory as shown below.
1. ./data/R1A/R1A.prepin
2. ./data/R1A/frcmod.R1A

Next, we use the generated parameters for 3D structure modeling, feature extraction, machine learning, and other operations.

1. Generate the 3D structure of the glycosylated arginine peptide
```python
from stapep.structure import Structure

# Import the additional_residues to specify the path of the parameter file
# R1A is the name of the non-standard amino acid
# data/R1A/R1A.prepin is the path of the prepin file
# data/R1A/frcmod.R1A is the path of the frcmod file
additional_residues = {
    'R1A': (
        'data/R1A/R1A.prepin',
        'data/R1A/frcmod.R1A',
    )
}

seq = 'Ac-BATPRRR-BLBR-R1A-FKRLQ' # Define the peptide sequence
st = Structure(verbose=True) # Initialize the structure class

st.de_novo_3d_structure(seq=seq, 
                        output_pdb='data/R1A_peptide.pdb', 
                        additional_residues=additional_residues)
```
![R1A_peptide](https://github.com/dahuilangda/stapep_package/blob/master/stapep/example/img/peptide_within_non_standard_aa.png)

2. Molecular dynamics simulation of the glycosylated arginine peptide
```python
from stapep.molecular_dynamics import PrepareProt, Simulation

output = 'data/R1A'
additional_residues = {
    'R1A': (
        'data/R1A/R1A.prepin',
        'data/R1A/frcmod.R1A',
    )
}

# Prepare the protein for molecular dynamics simulation
pp = PrepareProt(seq, output, method='alphafold', additional_residues=additional_residues)
pp._gen_prmtop_and_inpcrd_file()

sim = Simulation(output)
sim.setup(
        type='implicit', # 'explicit' or 'implicit'
        solvent='water', # 'water' or 'chloroform'
        temperature=300, # Kelvin
        friction=1, # ps^-1
        timestep=2, # fs
        interval=100, # ps
        nsteps=500000 # 1 ns
    )
sim.minimize()
sim.run()
```

3. Feature extraction of trajectory of molecular dynamics simulation of the glycosylated arginine peptide
```python
import os
from stapep.utils import PhysicochemicalPredictor

# initialize the PhysicochemicalPredictor class, which automatically loads trajectory using pytraj.
pcp = PhysicochemicalPredictor(sequence=seq, 
                                topology_file=os.path.join(output, 'pep_vac.prmtop'),
                                trajectory_file=os.path.join(output, 'traj.dcd'),
                                start_frame=0)

# You can call the calc_helix_percent, calc_extend_percent, and calc_loop_percent methods 
# to respectively calculate the percentages of helix, strand, and coil in the secondary structure.

print('helix percent: ', pcp.calc_helix_percent())
print('sheet percent: ', pcp.calc_extend_percent())
print('loop percent: ', pcp.calc_loop_percent())

# save the mean structure of the trajectory
pcp._save_mean_structure(os.path.join(output, 'mean_structure.pdb'))

# calculate the weight of the protein using the mean structure
print('weight', pcp.calc_weight(os.path.join(output, 'mean_structure.pdb')))
# calculate the Mean B-factor, Molecular Surface, Mean Gyration Radius, Hydrophobic Index, and 3D-PSA
print('mean bfactor: ', pcp.calc_mean_bfactor())
print('mol surf: ', pcp.calc_mean_molsurf())
print('mean gyrate: ', pcp.calc_mean_gyrate())
print('hydrophobic index: ', pcp.calc_hydrophobic_index(os.path.join(output, 'mean_structure.pdb')))
print('psa: ', pcp.calc_psa(os.path.join(output, 'mean_structure.pdb')))
print('total number of hydrogen bonds: ', pcp.calc_n_hbonds())

# Output:
helix percent:  0.20142222222222222
sheet percent:  4.444444444444444e-05
loop percent:  0.7985333333333332
weight 2414.0389999999984
mean bfactor:  190.2779644700375
mol surf:  1830.2296376834504
mean gyrate:  9.026684269192096
hydrophobic index:  -0.7733333333333332
psa:  1002.5029907226562
total number of hydrogen bonds:  9
```

4. Feature extraction of sequence of the glycosylated arginine peptide
```python
from stapep.utils import ProtParamsSeq
# initialize the ProtParamsSeq class
pps = ProtParamsSeq(seq)
# You can call the following methods to calculate the physicochemical properties of the peptide.
# length, weight, hydrophobicity index, charge, charge density, aromaticity, fraction of arginine, fraction of lysine, lyticity index, and isoelectric point.
print('length: ', pps.seq_length)
# print('weight: ', pps.weight) # This method is not implemented yet if the sequence contains non-standard amino acids.
# print('hydrophobicity index: ', pps.hydrophobicity_index) # This method is not implemented yet if the sequence contains non-standard amino acids.
print('charge', pps.calc_charge(pH=7.0)) # Be careful if the sequence contains non-standard amino acids, the charge will be calculated based on the standard amino acids.
# print('charge_density', pps.calc_charge_density(pH=7.0)) # The charge density is calculated by dividing the charge by the length of the peptide.
print('aromaticity', pps.aromaticity) # Be careful if the sequence contains non-standard amino acids, the aromaticity will be calculated based on the standard amino acids.
print('fraction_arginine', pps.fraction_arginine)
print('fraction_lysine', pps.fraction_lysine)
# print('lyticity index: ', pps.lyticity_index) # This method is not implemented yet if the sequence contains non-standard amino acids.
print('isoelectric_point: ', pps.isoelectric_point) # Be careful if the sequence contains non-standard amino acids, the isoelectric point will be calculated based on the standard amino acids.

# Output:
length:  19
charge 6.996
aromaticity 0.05263157894736842
fraction_arginine 0.3157894736842105
fraction_lysine 0.05263157894736842
isoelectric_point:  11.999967765808105
```
_Demo: examples/ex4_non_standard_amino_acid.ipynb_