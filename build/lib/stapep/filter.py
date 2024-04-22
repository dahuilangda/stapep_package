import os
import shutil
from collections import defaultdict
from itertools import product
import pandas as pd 

from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
from Bio import SeqIO

from stapep.utils import SeqPreProcessing
from stapep.params import amino_acid_dict


class Filter(object):
    '''
        Filter sequences based on some criteria.

        Reference: Söylemez, Ümmü Gülsüm, et al. 
                  "Prediction of Linear Cationic Antimicrobial Peptides Active against Gram-Negative and Gram-Positive Bacteria Based on Machine Learning Models." 
                  Applied Sciences 12.7 (2022): 3631.

        Params:
            df: pd.DataFrame
                DataFrame of sequences
            name_col: str
                Name of column that contains sequence
            seq_col: str
                Column name of sequence
    '''

    def __init__(self, 
                 df: pd.DataFrame, 
                 name_col: str='name', 
                 seq_col: str='sequence'):
        self.df = df
        self.name_col = name_col
        self.seq_col = seq_col
        self.seqpp = SeqPreProcessing()

    def _check_col_exist(self):
        '''
            Check if columns exist
        '''
        if self.name_col not in self.df.columns:
            raise ValueError(f'{self.name_col} is not a valid column name')
        if self.seq_col not in self.df.columns:
            raise ValueError(f'{self.seq_col} is not a valid column name')

    def _filter_seq_length(self, seq: str, min_length: int=10, max_length: int=50) -> str:
        '''
            Filter sequence by length
            return: bool (True/False)
        '''
        seq_list = self.seqpp._seq_to_list(seq=seq)
        return len(seq_list) >= min_length and len(seq_list) <= max_length

    def filter_seq_length(self, min_length: int=10, max_length: int=50) -> pd.DataFrame:
        '''
            Filter sequence by length
            return: pd.DataFrame
        '''
        self.df['seq_length_valid'] = self.df[self.seq_col].apply(self._filter_seq_length, args=(min_length, max_length,))
        return self.df[self.df['seq_length_valid'] == True]

    def _df_to_fasta(self, output: str) -> str:
        '''
            Convert seq of df to fasta format

            Params:
                output: output directory
        '''
        fasta = [f'>{name}\n{seq}\n' for name, seq in zip(self.df[self.name_col], self.df[self.seq_col])]
        with open(f'{output}/seqs.fasta', 'w') as f:
            f.writelines(fasta)

        return f'{output}/seqs.fasta'

    def _fasta_to_df(self, output_file: str) -> pd.DataFrame:
        '''
            Convert fasta to df
            return: pd.DataFrame
        '''
        records = list(SeqIO.parse(output_file, 'fasta'))
        names = []
        seqs = []
        for record in records:
            names.append(str(record.id))
            seqs.append(str(record.seq))
        df = pd.DataFrame(columns=['name', 'sequence'])
        df['name'] = names
        df['sequence'] = seqs
        df = self._merge_score_to_new_df(df)
        return df

    def _merge_score_to_new_df(self, new_df: pd.DataFrame) -> pd.DataFrame:
        '''
            Merge score col of self.df to new_df according to seq col

            Args:
                new_df: pd.DataFrame
                    New DataFrame

            Return: pd.DataFrame
        '''

        new_df = pd.concat([new_df, self.df], axis=1, join='inner')
        if 'Unnamed: 0' in new_df.columns:
            del new_df['Unnamed: 0']
        return new_df

    def _check_cdhit_exist(self):
        '''
            Check if CD-HIT exist.
        '''
        if shutil.which('cd-hit') is None:
            raise FileNotFoundError('cd-hit is not installed; please install it first! sudo apt install cd-hit')

    def filter_identity_threshold(self, 
                                  threshold: float=0.7,
                                  temp_dir: str='/tmp') -> pd.DataFrame:
        '''
            Filter sequence by similarity threshold using CD-HIT program.

            Args:
                threshold: float
                    Threshold of similarity
                temp_dir: str
                    Temporary directory to store intermediate files

            Return: pd.DataFrame
        '''
        self._check_cdhit_exist()
        output = os.path.join(temp_dir, f'cdhit_{threshold}')
        if not os.path.exists(output):
            os.makedirs(output)
        fasta = self._df_to_fasta(str(output))
        fasta_out = f'{output}/seqs_filtered.fasta'

        if 0.7 < threshold <= 1.0:
            word_size = 5
        elif 0.6 < threshold <= 0.7:
            word_size = 4
        elif 0.5 < threshold <= 0.6:
            word_size = 3
        elif 0.4 < threshold <= 0.5:
            word_size = 2
        else:
            raise ValueError('threshold should be between 0.4 and 1.0')

        cmd = f'cd-hit -i {fasta} -o {fasta_out} -c {threshold} -n {word_size} -M 16000 -d 0 -T 8'
        # print(cmd)
        os.system(cmd)
        if os.path.exists(fasta_out):
            return  self._fasta_to_df(fasta_out)
        return pd.DataFrame()

    def drop_duplicate(self) -> pd.DataFrame:
        '''
            Drop duplicate sequences
            Return: pd.DataFrame
        '''
        seq_list = self.df[self.seq_col].tolist()
        seq_norm_list = [''.join(self.seqpp._seq_to_list(seq)) for seq in seq_list]
        self.df['seq_norm'] = seq_norm_list
        self.df.drop_duplicates(subset=['seq_norm'], inplace=True)
        del self.df['seq_norm']
        return self.df
    
    def _compute_kmer_frequencies(self, seq: str, k: int = 3) -> list:
        """
        Compute k-mer frequencies for the given sequence.

        Args:
            seq: str
                Sequence to compute k-mer frequencies
            k: int
                Length of k-mers

        Return: list
            List of k-mer frequencies
        """
        amino_acids = ''.join(amino_acid_dict.values())
        kmers = [''.join(p) for p in product(amino_acids, repeat=k)]
        kmer_counts = defaultdict(int)
        total_kmers = len(seq) - k + 1

        for i in range(total_kmers):
            kmer = seq[i:i+k]
            kmer_counts[kmer] += 1

        return [kmer_counts[kmer] / total_kmers for kmer in kmers]
    
    def cluster_sequences(self, cluster_col: str='cluster', n_clusters: int = 5, k: int = 3) -> pd.DataFrame:
        """
        Cluster sequences based on their k-mer frequencies.

        Args:
            cluster_col: str
                Name of the column containing the cluster labels
            n_clusters: int
                Number of clusters to form
            k: int
                Length of k-mers
        
        Return: pd.DataFrame
            DataFrame with an additional column 'cluster' containing the cluster labels
        """
        # Compute k-mer frequencies for each sequence
        self.df['kmer_frequencies'] = self.df[self.seq_col].apply(self._compute_kmer_frequencies, args=(k,))

        # Perform k-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=0)
        self.df[cluster_col] = kmeans.fit_predict(list(self.df['kmer_frequencies']))

        # Remove the 'kmer_frequencies' column
        self.df.drop(columns=['kmer_frequencies'], inplace=True)

        return self.df
    
    def split_dataset_by_cluster(self, df: pd.DataFrame, cluster_col: str = 'cluster', test_size: float = 0.2, random_state: int = 42):
        """
        Split dataset into train and test sets based on the clustering results.

        Params:
            df: pd.DataFrame
                DataFrame containing the sequences and their cluster labels
            cluster_col: str
                Name of the column containing the cluster labels
            test_size: float
                Proportion of the dataset to include in the test split
            random_state: int
                Random seed used for reproducibility

        Return:
            train_df: pd.DataFrame
                DataFrame containing the training set
            test_df: pd.DataFrame
                DataFrame containing the test set
        """
        train_df = pd.DataFrame(columns=df.columns)
        test_df = pd.DataFrame(columns=df.columns)

        for cluster_label in df[cluster_col].unique():
            cluster_data = df[df[cluster_col] == cluster_label]
            cluster_train, cluster_test = train_test_split(cluster_data, test_size=test_size, random_state=random_state)
            train_df = pd.concat([train_df, cluster_train], ignore_index=True)
            test_df = pd.concat([test_df, cluster_test], ignore_index=True)

        return train_df, test_df
    
if __name__ == '__main__':
    df = pd.read_csv('../example/datasets/Stapled-peptide_permeability.csv')
    
    filter = Filter(df, seq_col='SEQUENCE', name_col='name')
    print('Original', df.shape)
    df = filter.drop_duplicate()
    print('After drop duplicate', df.shape)

    filter = Filter(df, seq_col='SEQUENCE', name_col='name')
    df = filter.filter_seq_length()
    print('After filter seq length', df.shape)

    filter = Filter(df, seq_col='SEQUENCE', name_col='name')
    df = filter.filter_identity_threshold(threshold=0.95)
    print('After filter identity threshold', df.shape)

    filter = Filter(df, seq_col='SEQUENCE', name_col='name')
    df = filter.cluster_sequences(n_clusters=5, k=2)

    print(df.value_counts('cluster'))

    filter = Filter(df, seq_col='SEQUENCE', name_col='name')
    train_df, test_df = filter.split_dataset_by_cluster(df,
                                                        cluster_col='cluster',
                                                        test_size=0.2,
                                                        random_state=42)
    print('Train', train_df.shape)
    print('Test', test_df.shape)