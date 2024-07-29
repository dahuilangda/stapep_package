import shutil
import numpy as np
import pandas as pd
import os
import re

from functools import wraps

from Bio import SeqIO
import pandas as pd
import pytraj as pt
import signal, functools
import networkx as nx
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors

from stapep.params import amino_acid_dict, reversed_amino_acid_dict, li_dict, weight_dict, hydrophobic_dict, hydrophilic_residues

def proxy_decorator(proxy):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            os.environ['http_proxy'] = proxy
            os.environ['https_proxy'] = proxy
            result = func(*args, **kwargs)
            os.environ['http_proxy'] = ''
            os.environ['https_proxy'] = ''
            return result
        return wrapper
    return decorator

def timeout(seconds, error_message="Timeout Error: the cmd have not finished."):
    def decorated(func):
        result = 'Unknow'
        def _handle_timeout(signum, frame):
            global result
            result = error_message
            raise TimeoutError(error_message)
 
        def wrapper(*args, **kwargs):
            global result
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
 
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
                return result
            return result

        return functools.wraps(func)(wrapper)
    return decorated


class PermDataSet(object):
    '''
        Permutation dataset from hundreds of papers.
    '''

class ProtParamsSeq(object):
    '''
        Calculate some params based on sequence...

        attributes:
            1. seq_length: length of sequence
            2. weight: return the weight of sequence.
            3. hydrophobicity_index: return the hydrophobicity index of sequence.
            4. lyticity_index: return the lytic index of sequence.
            5. charge: return the charge of sequence.
            6. charge_density: return the charge density of sequence.
            7. fraction_arginine: return the fraction of arginine in sequence.
            8. fraction_lysine: return the fraction of lysine in sequence.
            9. isoelectric_point: return the isoelectric point of sequence.

        Methods:
            1. plot_lyticity_index: plot the lytic index of sequence, and return the lyticity index.
    '''
    def __init__(self, seq: str, additional_params=None) -> None:
        '''
            Args:
                seq: sequence of peptide.
                additional_params: additional residues to update li_dict and hydrophobic_dict.
                    eg: {
                            'lyticity_index': {'Aib': 8.493, 'NLE': 24.442},
                            'hydrophobicity_index': {'Aib': 8.493, 'NLE': 24.442},
                            'weight': {'Aib': 8.493, 'NLE': 24.442},
                            'hydrophilic_residues': {'Aib': False, 'NLE': True},
                        }
        '''
        
        additional_residues = {}
        if additional_params is not None:
            if 'lyticity_index' in additional_params.keys():
                for k, v in additional_params['lyticity_index'].items():
                    if v:
                        additional_residues[k] = None
            if 'hydrophobicity_index' in additional_params.keys():
                for k, v in additional_params['hydrophobicity_index'].items():
                    if v:
                        additional_residues[k] = None
            if 'weight' in additional_params.keys():
                for k, v in additional_params['weight'].items():
                    if v:
                        additional_residues[k] = None
            if 'hydrophilic_residues' in additional_params.keys():
                for k, v in additional_params['hydrophilic_residues'].items():
                    if v:
                        additional_residues[k] = None


        self.seq = seq
        seqpp = SeqPreProcessing(additional_residues=additional_residues)
        self.seq_to_list = seqpp._seq_to_list(self.seq)
        self.weight_dict = weight_dict
        self.li_dict = li_dict
        self.hydrophobic_dict = hydrophobic_dict
        self.hydrophilic_residues = hydrophilic_residues

        if additional_params is not None:
            if 'lyticity_index' in additional_params.keys():
                self.li_dict.update(additional_params['lyticity_index'])
            if 'hydrophobicity_index' in additional_params.keys():
                self.hydrophobic_dict.update(additional_params['hydrophobicity_index'])
            if 'weight' in additional_params.keys():
                self.weight_dict.update(additional_params['weight'])
            if 'hydrophilic_residues' in additional_params.keys():
                for k, v in additional_params['hydrophilic_residues'].items():
                    if v:
                        self.hydrophilic_residues.append(k)

    @property
    def seq_length(self) -> int:
        '''
            Get sequence length.
        '''
        seq_length = len(self.seq_to_list)
        if 'Ac' in self.seq_to_list:
            seq_length -= 1
        if 'NH2' in self.seq_to_list:
            seq_length -= 1
        return seq_length

    @property
    def hydrophobicity_index(self) -> float:
        '''
            Calculate Kyte & Doolittle index of hydrophobicity of sequence.
            DOI: 10.1016/0022-2836(82)90515-0.
        '''
        for aa in self.seq_to_list:
            if aa not in self.hydrophobic_dict.keys() and aa not in ['Ac', 'NH2']:
                raise ValueError(f'{aa} is not a valid amino acid.')

        hydrophobic_index_list = [self.hydrophobic_dict[aa] for aa in self.seq_to_list if aa not in ['Ac', 'NH2']]
        return np.mean(hydrophobic_index_list)

    @property
    def weight(self) -> float:
        '''
            Calculate the weight of peptide.
        '''
        for aa in self.seq_to_list:
            if aa not in weight_dict.keys() and aa not in ['Ac', 'NH2']:
                raise ValueError(f'{aa} is not a valid amino acid.')

        weights = weight_dict
        water = 18.02
        carbon = 12.01
        _ACE = 41.01
        _NME = 0.023
        num_stapled_aa = len([aa for aa in self.seq_to_list if aa in ['S3', 'S5', 'S8', 'R3', 'R5', 'R8']])
        _ac = 1 if 'Ac' in self.seq_to_list else 0
        _nh2 = 1 if 'NH2' in self.seq_to_list else 0
        weight_list = [weights[aa] for aa in self.seq_to_list if aa not in ['Ac', 'NH2']]
        return np.sum(weight_list) - (len(self.seq_to_list)-1) * water - num_stapled_aa * carbon + _ac*_ACE + _nh2*_NME

    @property
    def lyticity_index(self) -> float:
        '''
            Mourtada, Rida, et al. 
            "Design of stapled antimicrobial peptides that are stable, nontoxic and kill antibiotic-resistant bacteria in mice." 
            Nature biotechnology 37.10 (2019): 1186-1197.

            Return the lyticity index of sequence.

            copy from https://www.walenskylab.org/hnm-landing
        '''
        sequence = self.seq_to_list
        
        for aa in sequence:
            if aa not in self.li_dict.keys() and aa not in ['Ac', 'NH2']:
                raise ValueError(f'{aa} is not a valid amino acid.')

        sequence = [aa for aa in sequence if aa not in ['Ac', 'NH2']]

        # Assign lyticity values to each residue in the peptide sequence
        lyticity_assignment = [self.li_dict[aa] for aa in sequence]
        
        # Make array of i+4 sums
        i_plus_4_sums = []
        for i in range(len(sequence) - 4):
            if sequence[i] not in self.hydrophilic_residues and sequence[i + 4] not in self.hydrophilic_residues:
                i_plus_4_sums.append(lyticity_assignment[i] + lyticity_assignment[i + 4])
        
        # Make array of i+3 sums
        i_plus_3_sums = []
        for i in range(len(sequence) - 3):
            if sequence[i] not in self.hydrophilic_residues and sequence[i + 3] not in self.hydrophilic_residues:
                i_plus_3_sums.append(lyticity_assignment[i] + lyticity_assignment[i + 3])
        
        # Sum everything together to get the final lyticity value
        total_lyticity = sum(i_plus_4_sums) + sum(i_plus_3_sums)
        
        return total_lyticity

    # @property
    # def lyticity_index(self) -> float:
    #     '''
    #         Mourtada, Rida, et al. 
    #         "Design of stapled antimicrobial peptides that are stable, nontoxic and kill antibiotic-resistant bacteria in mice." 
    #         Nature biotechnology 37.10 (2019): 1186-1197.

    #         Return the lyticity index of sequence.

    #     '''
    #     for aa in self.seq_to_list:
    #         if aa not in self.li_dict.keys() and aa not in ['Ac', 'NH2']:
    #             raise ValueError(f'{aa} is not a valid amino acid.')

    #     width_list = []
    #     nodes = [s for s in self.seq_to_list if s in self.li_dict]
    #     h_dict = dict(self.li_dict)
    #     for idx in range(len(nodes)-3):
    #         if (nodes[idx].split('_')[0] in self.li_dict) and (nodes[idx+3].split('_')[0] in self.li_dict):
    #             width_list.append((h_dict[nodes[idx].split('_')[0]]+h_dict[nodes[idx+3].split('_')[0]])/5)
    #         if idx + 4 < len(nodes) and (nodes[idx].split('_')[0] in self.li_dict) and (nodes[idx + 4].split('_')[0] in self.li_dict):
    #             width_list.append((h_dict[nodes[idx].split('_')[0]]+h_dict[nodes[idx+4].split('_')[0]])/5)
    #     return np.sum(width_list)

    def calc_charge(self, pH: float=7.0, amide: bool=False) -> float:
        # https://github.com/alexarnimueller/modlAMP/blob/master/modlamp/descriptors.py
        """Calculates charge of a single sequence. The method used is first described by Bjellqvist. In the case of
        amidation, the value for the  'Cterm' pKa is 15 (and Cterm is added to the pos_pks dictionary.
        The pKa scale is extracted from: http://www.hbcpnetbase.com/ (CRC Handbook of Chemistry and Physics, 96th ed).

        **pos_pks** = {'Nterm': 9.38, 'K': 10.67, 'R': 12.10, 'H': 6.04}

        **neg_pks** = {'Cterm': 2.15, 'D': 3.71, 'E': 4.15, 'C': 8.14, 'Y': 10.10}

        args:
            pH: {float} pH at which to calculate peptide charge.
            amide: {boolean} whether the sequences have an amidated C-terminus.

        return: {array} descriptor values in the attribute :py:attr:`descriptor
        """
        
        def _count_aas() -> dict:
            _aa_count_set = ['K', 'R', 'H', 'D', 'E', 'C', 'Y']
            return {aa: self.seq_to_list.count(aa) for aa in _aa_count_set}

        if amide:
            pos_pks = {'Nterm': 9.38, 'K': 10.67, 'R': 12.10, 'H': 6.04}
            neg_pks = {'Cterm': 15., 'D': 3.71, 'E': 4.15, 'C': 8.14, 'Y': 10.10}
        else:
            pos_pks = {'Nterm': 9.38, 'K': 10.67, 'R': 12.10, 'H': 6.04}
            neg_pks = {'Cterm': 2.15, 'D': 3.71, 'E': 4.15, 'C': 8.14, 'Y': 10.10}
        
        aa_content = _count_aas()
        aa_content['Nterm'] = 1.0
        aa_content['Cterm'] = 1.0
        pos_charge = 0.0
        for aa, pK in pos_pks.items():
            c_r = 10 ** (pK - pH)
            partial_charge = c_r / (c_r + 1.0)
            pos_charge += aa_content[aa] * partial_charge
        neg_charge = 0.0
        for aa, pK in neg_pks.items():
            c_r = 10 ** (pH - pK)
            partial_charge = c_r / (c_r + 1.0)
            neg_charge += aa_content[aa] * partial_charge
        return round(pos_charge - neg_charge, 3)
    
    def calc_charge_density(self, pH: float=7.0, amide: bool=False) -> float:
        '''
            Calculate the charge density of sequence. Which is the charge of sequence divided by its weight.

            args:
                pH: pH value
                amide: whether the sequences have an amidated C-terminus.
            return:
                charge density
        '''
        return self.calc_charge(pH, amide) / self.weight

    @property
    def fraction_arginine(self) -> float:
        '''
            Return the fraction of arginine in sequence.
        '''
        num_arg = len([aa for aa in self.seq_to_list if aa == 'R'])
        return num_arg / self.seq_length

    @property
    def fraction_lysine(self) -> float:
        '''
            Return the fraction of lysine in sequence.
        '''
        num_arg = len([aa for aa in self.seq_to_list if aa == 'K'])
        return num_arg / self.seq_length

    @property
    def aromaticity(self):
        '''
            Return the aromaticity of sequence.
        '''
        def _count_aas() -> dict:
            _aa_count_set = ['F', 'W', 'Y']
            return {aa: self.seq_to_list.count(aa) for aa in _aa_count_set}

        aa_content = _count_aas()
        return (aa_content['F'] + aa_content['W'] + aa_content['Y']) / self.seq_length

    def plot_lyticity_index(self, output_path: str) -> float:
        '''
            Mourtada, Rida, et al. 
            "Design of stapled antimicrobial peptides that are stable, nontoxic and kill antibiotic-resistant bacteria in mice." 
            Nature biotechnology 37.10 (2019): 1186-1197.

            Args:
                output_path: save path of plot.
            
            Return:
                lyticity_index: lyticity index of sequence.
        '''

        sequence = [aa for aa in self.seq_to_list if aa not in ['Ac', 'NH2']]

        G = nx.Graph()
        h_dict = dict(self.li_dict)
        nodes = [f'{s}_{step}' for step,s in enumerate(sequence) if s in self.li_dict]
        labels = {f'{s}_{step}':r'$%s_{%s}$'%(s, step+1) for step, s in enumerate(sequence) if s in self.li_dict}
        
        connected_nodes = set()
        width_list = []
        
        # Create edges and calculate widths
        for idx in range(len(nodes) - 3):
            node1 = nodes[idx].split('_')[0]
            node2 = nodes[idx+3].split('_')[0]
            if node1 in self.li_dict and node2 in self.li_dict and node1 not in self.hydrophilic_residues and node2 not in self.hydrophilic_residues:
                G.add_edge(nodes[idx], nodes[idx+3])
                width_list.append((h_dict[node1] + h_dict[node2]) / 5)
                connected_nodes.add(nodes[idx])
                connected_nodes.add(nodes[idx+3])
                
            if idx + 4 < len(nodes):
                node3 = nodes[idx+4].split('_')[0]
                if node1 in self.li_dict and node3 in self.li_dict and node1 not in self.hydrophilic_residues and node3 not in self.hydrophilic_residues:
                    G.add_edge(nodes[idx], nodes[idx+4])
                    width_list.append((h_dict[node1] + h_dict[node3]) / 5)
                    connected_nodes.add(nodes[idx])
                    connected_nodes.add(nodes[idx+4])
        
        # Only add nodes that have connections
        G = G.subgraph(connected_nodes)
        
        position = nx.spring_layout(G)
        fig, ax = plt.subplots(1, 1, figsize=(20, 12))
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        nx.draw_networkx_nodes(G, position, node_color='black', alpha=1.0, ax=ax, node_size=2200)
        nx.draw_networkx_nodes(G, position, node_color='bisque', alpha=1.0, ax=ax, node_size=2000)
        nx.draw_networkx_edges(G, position, edge_color='orangered', alpha=1, ax=ax, width=width_list)
        nx.draw_networkx_labels(G, position, {node: labels[node] for node in G.nodes()}, font_size=16, ax=ax)
        
        li = sum(width_list) * 5
        # li = self.lyticity_index
        plt.text(0.1, 1.0, 'Lyticity index: ' + '%.2f' % li, fontsize=20, ha='center', va='center', transform=ax.transAxes)
        plt.savefig(f'{output_path}', format='SVG')
        return li


    def plot_lyticity_index_interactive(self, output_path: str) -> float:

        '''
            Mourtada, Rida, et al. 
            "Design of stapled antimicrobial peptides that are stable, nontoxic and kill antibiotic-resistant bacteria in mice." 
            Nature biotechnology 37.10 (2019): 1186-1197.

            Args:
                output_path: save path of plot.
            
            Return:
                lyticity_index: lyticity index of sequence.
        '''
    
        try:
            import networkx as nx
            from pyvis.network import Network
        except ImportError:
            raise ImportError('Please install networkx and pyvis to use this function. eg: pip install networkx pyvis')

        sequence = [aa for aa in self.seq_to_list if aa not in ['Ac', 'NH2']]

        G = nx.Graph()
        h_dict = dict(self.li_dict)
        nodes = [f'{s}_{step}' for step, s in enumerate(sequence) if s in self.li_dict]
        labels = {f'{s}_{step}': f'{s}{step+1}' for step, s in enumerate(sequence) if s in self.li_dict}
        
        connected_nodes = set()
        width_list = []
        
        # Create edges and calculate widths
        for idx in range(len(nodes) - 3):
            node1 = nodes[idx].split('_')[0]
            node2 = nodes[idx+3].split('_')[0]
            if node1 in self.li_dict and node2 in self.li_dict and node1 not in self.hydrophilic_residues and node2 not in self.hydrophilic_residues:
                G.add_edge(nodes[idx], nodes[idx+3])
                width_list.append((h_dict[node1] + h_dict[node2]) / 5)
                connected_nodes.add(nodes[idx])
                connected_nodes.add(nodes[idx+3])
                
            if idx + 4 < len(nodes):
                node3 = nodes[idx+4].split('_')[0]
                if node1 in self.li_dict and node3 in self.li_dict and node1 not in self.hydrophilic_residues and node3 not in self.hydrophilic_residues:
                    G.add_edge(nodes[idx], nodes[idx+4])
                    width_list.append((h_dict[node1] + h_dict[node3]) / 5)
                    connected_nodes.add(nodes[idx])
                    connected_nodes.add(nodes[idx+4])
        
        # Only add nodes that have connections
        G = G.subgraph(connected_nodes)
        
        net = Network(height='750px', width='100%', notebook=True)
        
        # Add nodes with labels and sizes
        for node in G.nodes:
            net.add_node(node, label=labels[node], size=20, title=labels[node])
        
        # Add edges with widths
        for edge, width in zip(G.edges, width_list):
            net.add_edge(edge[0], edge[1], width=width)
        
        # Save the network visualization to an HTML file
        net.save_graph(output_path)

        li = sum(width_list) * 5
        return li


    # def plot_lyticity_index_3d(self, output_path: str) -> float:
    #     try:
    #         from plotly.graph_objects import Scatter3d, Layout, Figure
    #     except ImportError:
    #         raise ImportError('Please install plotly to use this function. eg: pip install plotly')

    #     sequence = [aa for aa in self.seq_to_list if aa not in ['Ac', 'NH2']]
    #     G = nx.Graph()
    #     h_dict = dict(self.li_dict)
    #     nodes = [f'{s}_{step}' for step, s in enumerate(sequence) if s in self.li_dict]
        
    #     # Use HTML-like syntax for subscripts in labels
    #     labels = {f'{s}_{step}': f'{s}<sub>{step + 1}</sub>' for step, s in enumerate(sequence) if s in self.li_dict}
        
    #     connected_nodes = set()
    #     width_list = []
        
    #     # Create edges and calculate widths
    #     for idx in range(len(nodes) - 3):
    #         node1 = nodes[idx].split('_')[0]
    #         node2 = nodes[idx + 3].split('_')[0]
    #         if node1 in self.li_dict and node2 in self.li_dict and node1 not in self.hydrophilic_residues and node2 not in self.hydrophilic_residues:
    #             G.add_edge(nodes[idx], nodes[idx + 3])
    #             width_list.append((h_dict[node1] + h_dict[node2]) / 5)
    #             connected_nodes.add(nodes[idx])
    #             connected_nodes.add(nodes[idx + 3])
                
    #         if idx + 4 < len(nodes):
    #             node3 = nodes[idx + 4].split('_')[0]
    #             if node1 in self.li_dict and node3 in self.li_dict and node1 not in self.hydrophilic_residues and node3 not in self.hydrophilic_residues:
    #                 G.add_edge(nodes[idx], nodes[idx + 4])
    #                 width_list.append((h_dict[node1] + h_dict[node3]) / 5)
    #                 connected_nodes.add(nodes[idx])
    #                 connected_nodes.add(nodes[idx + 4])
        
    #     # Only add nodes that have connections
    #     G = G.subgraph(connected_nodes)
        
    #     position = nx.spring_layout(G, dim=3)
        
    #     edge_traces = []
    #     for i, edge in enumerate(G.edges()):
    #         x0, y0, z0 = position[edge[0]]
    #         x1, y1, z1 = position[edge[1]]
    #         edge_traces.append(
    #             Scatter3d(
    #                 x=[x0, x1, None],
    #                 y=[y0, y1, None],
    #                 z=[z0, z1, None],
    #                 line=dict(width=width_list[i], color='orangered'),
    #                 mode='lines'
    #             )
    #         )
        
    #     node_trace = Scatter3d(
    #         x=[position[node][0] for node in G.nodes()],
    #         y=[position[node][1] for node in G.nodes()],
    #         z=[position[node][2] for node in G.nodes()],
    #         mode='markers+text',
    #         text=[labels[node] for node in G.nodes()],
    #         hoverinfo='text',
    #         marker=dict(
    #             showscale=False,
    #             color='bisque',
    #             size=10,
    #             line=dict(width=2, color='black')),
    #         textposition='middle center'
    #     )
        
    #     layout = Layout(
    #         showlegend=False,
    #         margin=dict(b=0, l=0, r=0, t=0),
    #         scene=dict(
    #             xaxis=dict(showbackground=False, showticklabels=False),
    #             yaxis=dict(showbackground=False, showticklabels=False),
    #             zaxis=dict(showbackground=False, showticklabels=False)
    #         ),
    #         annotations=[dict(
    #             text=f'Lyticity index: {sum(width_list) * 5:.2f}',
    #             showarrow=False,
    #             xref="paper", yref="paper",
    #             x=0.5, y=1.05,
    #             font=dict(size=16)
    #         )]
    #     )
        
    #     fig = Figure(data=[node_trace] + edge_traces, layout=layout)
    #     fig.write_html(output_path)
    #     return sum(width_list) * 5



    # def plot_lyticity_index(self, output_path: str) -> float:
    #     '''
    #         Mourtada, Rida, et al. 
    #         "Design of stapled antimicrobial peptides that are stable, nontoxic and kill antibiotic-resistant bacteria in mice." 
    #         Nature biotechnology 37.10 (2019): 1186-1197.

    #         Args:
    #             output_path: save path of plot.
            
    #         Return:
    #             lyticity_index: lyticity index of sequence.
    #     '''

    #     G = nx.Graph()
    #     h_dict = dict(self.li_dict)
    #     nodes = [f'{s}_{step}' for step,s in enumerate(self.seq_to_list) if s in self.li_dict]
    #     labels = {f'{s}_{step}':r'$%s_{%s}$'%(s, step+1) for step, s in enumerate(self.seq_to_list) if s in self.li_dict}
    #     G.add_nodes_from(nodes)
    #     width_list = []
    #     for idx in range(len(nodes)-3):
    #         if (nodes[idx].split('_')[0] in self.li_dict) and (nodes[idx+3].split('_')[0] in self.li_dict):
    #             print(f'{nodes[idx]} {nodes[idx+3]}')
    #             G.add_edge(nodes[idx], nodes[idx+3])
    #             width_list.append((h_dict[nodes[idx].split('_')[0]]+h_dict[nodes[idx+3].split('_')[0]])/5)
    #         if idx + 4 < len(nodes) and (nodes[idx].split('_')[0] in self.li_dict) and (nodes[idx + 4].split('_')[0] in self.li_dict):
    #             print(f'{nodes[idx]} {nodes[idx+4]}')
    #             G.add_edge(nodes[idx], nodes[idx+4])
    #             width_list.append((h_dict[nodes[idx].split('_')[0]]+h_dict[nodes[idx+4].split('_')[0]])/5)

    #     position = nx.spring_layout(G)
    #     fig, ax = plt.subplots(1, 1, figsize=(20, 12))
    #     ax.spines['right'].set_visible(False)
    #     ax.spines['left'].set_visible(False)
    #     ax.spines['top'].set_visible(False)
    #     ax.spines['bottom'].set_visible(False)
    #     nx.draw_networkx_nodes(G, position, node_color='black', alpha = 1.0, ax=ax, node_size=2200)
    #     nx.draw_networkx_nodes(G, position, node_color='bisque', alpha = 1.0, ax=ax, node_size=2000)
    #     nx.draw_networkx_edges(G, position, edge_color='orangered', alpha = 1, ax=ax, width=width_list)
    #     nx.draw_networkx_labels(G, position, labels, font_size=16, ax=ax)
    #     li = sum(width_list*5)
    #     plt.text(0.1, 1.0, 'Lyticity index: '+'%.2f'%li, fontsize=20, ha='center', va='center', transform=ax.transAxes)
    #     plt.savefig(f'{output_path}', format='SVG')
    #     return li

    @property
    def isoelectric_point(self) -> float:
        '''
            Calculate the isoelectric point of sequence.

            Return:
                isoelectric_point: isoelectric point of sequence.
        '''
        seq_fixed = [x if len(x) == 1 else 'X' for x in self.seq_to_list ]
        seq_fixed = ''.join(seq_fixed)
        from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
        protein = IP(seq_fixed)
        return protein.pi()

class MDAnalysisHandler(object):
    '''
        MDAnalysis handler, Inherited by other classes.
    '''
    def __init__(self, 
                 topology_file: str, 
                 trajectory_file: str,
                 start_frame: int=500,
                 reimage: bool=False) -> None:
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.start_frame = start_frame
        self.reimage = reimage

    def _check_input_validity(self) -> None:
        '''
            Check input validity
        '''
        _, toplogy_ext = os.path.splitext(self.topology_file)
        if toplogy_ext not in ['.prmtop', '.top', '.pdb']:
            raise ValueError(f'{self.topology_file} is not a topology file. Expected .prmtop or .top or .pdb')

        _, trajectory_ext = os.path.splitext(self.trajectory_file)
        if trajectory_ext not in ['.dcd', '.xtc', '.trr', '.nc', '.mdcrd', '.binpos']:
            raise ValueError(f'{self.trajectory_file} is not a trajectory file. Expected .dcd or .xtc or .trr or .nc or .mdcrd or .binpos')

    def get_trajectory(self) -> pt.Trajectory:
        '''
            Load topology and trajectory
            Return pt.Trajectory
        '''
        self._check_input_validity()
        trajectory = pt.load(self.trajectory_file, top=self.topology_file)
        trajectory = trajectory[self.start_frame:]
        print(self.trajectory_file, self.topology_file)
        if self.reimage:
            pt.autoimage(trajectory) # Maybe very slowly!!!
        avg = pt.mean_structure(trajectory)
        pt.superpose(trajectory, ref=avg)
        return trajectory

    def save_traj(self, output_file: str, start_frame: int=500, reimage: bool=False) -> None:
        '''
            Save trajectory
        '''
        self.get_trajectory(start_frame=start_frame, reimage=reimage).save(output_file)


class PhysicochemicalPredictor(MDAnalysisHandler):
    '''
        Predict physicochemical properties of protein from MD trajectory:
        eg: b-factor, hydrophobicity, PSA, secondary structure, gyrate, etc.
    
    '''
    def __init__(self, sequence: str, 
                 topology_file: str, 
                 trajectory_file: str,
                 start_frame: int=500,
                 reimage: bool=False) -> None:
        super(PhysicochemicalPredictor, self).__init__(topology_file, trajectory_file, start_frame, reimage)
        self.sequence = sequence
        self.trajectory = self.get_trajectory()
    
    def _save_mean_structure(self, output_file: str) -> str:
        '''
            Save minimum RMSD structure along with average structure.
            TODO: how to calculate mean structure of MD trajectory? 
                  The minimum RMSD structure is the structure that is most similar to all the other structures in the trajectory, 
                  but it may not be representative of the entire trajectory.
        '''
        avg = pt.mean_structure(self.trajectory)
        rmsd_list = pt.rmsd(self.trajectory, ref=avg, mask='@CA')
        frame = self.trajectory[np.argmin(rmsd_list)]
        coord = frame.to_ndarray()
        coord = np.expand_dims(coord, axis=0)
        traj = pt.Trajectory(xyz=coord, top=self.trajectory.topology)
        traj.save(output_file)
        return output_file

    def extract_2d_structure(self, input_file:str) -> str:
        pep_ = Chem.MolFromPDBFile(input_file, removeHs=False)
        return Chem.MolToSmiles(pep_)

    def _get_bfactor(self) -> pd.DataFrame:
        '''
            Calculate b-factor from MD trajectory
            return: pd.DataFrame
                a numpy array of b-factor
                shape = (n_residues, 2)
        '''
        bfactor = pt.bfactors(self.trajectory, byres=True, mask='@CA')
        bfactor = pd.DataFrame(bfactor, columns=['residue', 'bfactor'])
        return bfactor

    def _calc_asa_for_average_structure(self, input_file: str) -> list[float]:
        '''
            Calculate relative ASA from average structure
            Return: list
        '''
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        self._save_mean_structure(input_file)

        p = PDBParser()
        structure = p.get_structure('PEPT', input_file)
        dssp = DSSP(structure[0], input_file)
        a_keys = list(dssp.keys())
        asa_list = []
        
        for a_key in a_keys:
            if dssp[a_key][3] == 'NA':
                asa_list.append(1.0) # stapled amino acid
            else:
                asa_list.append(float(dssp[a_key][3]))
        return asa_list

    def calc_hydrophobic_index(self, input_file: str, min_asa: float=0.2) -> list[float]:
        '''
            Calculate Kyte & Doolittle index of hydrophobicity of relative ASA > min_asa for average structure
            DOI: 10.1016/0022-2836(82)90515-0
        '''

        asa_list = self._calc_asa_for_average_structure(input_file)
        seqpp = SeqPreProcessing()

        seq_to_list = seqpp._seq_to_list(self.sequence)
        # if len(seq_to_list) != len(asa_list): # TODO: maybe is not correct
        #     raise ValueError(f'Sequence and ASA list length are not equal, {len(seq_to_list)} != {len(asa_list)}')

        hydrophobic_index_list = [hydrophobic_dict[aa] for aa, asa in zip(seq_to_list, asa_list) if asa > min_asa and aa in hydrophobic_dict]
        return np.mean(hydrophobic_index_list)

    def calc_n_hbonds(self) -> int:
        '''
            Calculate number of hydrogen bonds
        '''
        data = pt.search_hbonds(self.trajectory)
        return data.total_solute_hbonds()[0]

    def calc_mean_bfactor(self) -> float:
        '''
            Calculate mean b-factor
            return: float
        '''
        return self._get_bfactor().bfactor.mean()

    def _get_secondary_structure_matrix(self) -> pd.DataFrame:
        '''
            Calculate secondary structure from MD trajectory
            return: pd.DataFrame
        '''
        _, ss, _  = pt.dssp(self.trajectory, simplified=True)
        return pd.DataFrame(ss)

    def _calc_ss_percent(self, arg0):
        ss_matrix = self._get_secondary_structure_matrix()
        helix_list = []
        for idx in range(ss_matrix.shape[0]):
            count = len(ss_matrix.iloc[idx][ss_matrix.iloc[idx] == arg0])
            percent = count / float(ss_matrix.shape[1])
            helix_list.append(percent)
        return np.mean(helix_list)

    def calc_helix_percent(self) -> float:
        '''
            Calculate helix rate from secondary structure matrix
            return: float
        '''
        return self._calc_ss_percent('H')

    def calc_loop_percent(self) -> float:
        '''
            Calculate loop rate from secondary structure matrix
            return: float
        '''
        return self._calc_ss_percent('C')

    def calc_extend_percent(self) -> float:
        '''
            Calculate extend rate from secondary structure matrix
            return: float
        '''
        return self._calc_ss_percent('E')

    def _get_gyrate(self) -> list[float]:
        '''
            Calculate hydrophobicity from MD trajectory
            return: pd.DataFrame
        '''
        return pt.radgyr(self.trajectory)

    def calc_mean_gyrate(self):
        '''
            Calculate gyrate from MD trajectory
            return: float
        '''
        return np.mean(self._get_gyrate())

    def _get_molsurf(self) -> list[float]:
        '''
            Calculate molsurf from MD trajectory
            return: pd.DataFrame
        '''
        return pt.molsurf(self.trajectory)
    
    def calc_mean_molsurf(self) -> float:
        '''
            Calculate mean molsurf
            return: float
        '''
        return np.mean(self._get_molsurf()) # slowly method

    def calc_weight(self, input_file: str) -> float:
        '''
            Calculate weight from MD trajectory
            return: float
        '''
        return Descriptors.MolWt(Chem.MolFromSmiles(self.extract_2d_structure(input_file)))

    def calc_psa(self, input_file, obj_list=None, include_SandP = None, cmpd_name = None, atom_to_remove = None, **kwargs):    #TODO check that it is always consistently called 'psa3d' not 3dpsa or others
        """ ###TODO modify documentation
        Calculates the 3d polar surface area (3D-PSA) of molecules in Interface_Pymol for all the snapshots in a MD trajectory.
        (Contribution by Benjamin Schroeder)

        Parameters:
        -------------
        traj_file: str
            trajectory filename
        gro_file: str 
            coordinates filename (.gro or .pdb)
        cmpd_name: str, optional
            Name of the compound used as object name (default is "cmpd1")
        include_SandP: bool, optional 
            Set to False to exclude the S and P atoms from the calculation of the 3D-PSA. (Default = True) #TODO now you have S always included right?
        atom_to_remove: str, optional
            Single atom name of the atom to remove from the selection (Default = None). 
            Useful if you want to include only S or only P in the calculation of the 3D-PSA.
        Returns
        ----------
        dict_psa3d: dict
            Keys are mean (3d_psa_av), standard deviation (3d_psa_sd), and median (3d_psa_med) of the 3D-PSA calculated over the simulation time. 
            If cmpd_name is specified, it is returned in the dictionary.
        """
        try:
            from pymol import cmd
        except ImportError:
            raise ImportError('Extract 3D PSA is not possible beacause PyMol python handle not properly installed.')

        if cmpd_name is None:
            cmpd_name = "cmpd1" #TODO really needed?


        # Load trajectory and remove solvent and salts
        obj1 = cmpd_name
        cmd.reinitialize()
        cmd.load(input_file, object=obj1)
        # cmd.load_traj(traj_filename, object=obj1)

        atom_names = []
        cmd.iterate_state(-1, selection=f"{obj1} and not elem H", expression="atom_names.append(name)", space=locals())


        # IO
        obj_list = cmd.get_names("objects")
        assert len(obj_list) == 1
        # Select first (and only) object, which contains the trajectory for the solute molecule
        obj = obj_list[0]
        cmd.frame(0)
        states = range(1, cmd.count_states(obj) + 1)  # get all states of the object

        if include_SandP:
            select_string = f"(elem N or elem O or elem S or elem P or (elem H and (neighbor elem N+O+S+P))) and {obj} and not name {atom_to_remove}"

        else:
            select_string = f"(elem N or elem O or (elem H and (neighbor elem N+O))) and {obj} and not name {atom_to_remove}"

        # else:
        #     if include_SandP:
        #         select_string = "resn {} and (elem N or elem O or elem S or elem P or (elem H and (neighbor elem N+O+S+P))) and ".format(solute_resname) + obj   #@carmen add: "or elem S"
        #     else:
        #         select_string = "resn {} and (elem N or elem O or (elem H and (neighbor elem N+O))) and ".format(solute_resname) + obj  #@carmen add: "or elem S"
        ##Loop over all states
        psa = []
        for state in states:
            cmd.select("noh", select_string) #TODO is this really always called 'noh'?
            ###calc surface area
            # psa.append(float(cmd.get_area("noh", state=state)) * 0.01) #unit conversion from A^2 to nm^2
            psa.append(float(cmd.get_area("noh", state=state)))

        # return {f"ex4/{optimized}_psa3d" :  psa}
        return psa[0]


    def get_properties(self) -> pd.DataFrame:
        '''
            Get physicochemical properties from MD trajectory
            return: pd.DataFrame
        '''
        pass

class SeqPreProcessing(object):
    def __init__(self, additional_residues: dict=None):
        self.aa_dict = amino_acid_dict
        self.aa_reversed_dict = reversed_amino_acid_dict
        self.additional_residues = additional_residues

        if self.additional_residues is not None:
            if not isinstance(self.additional_residues, dict):
                raise ValueError('additional_residues must be a dict')

            # self.pattern_str = '(S3|S5|S8|R3|R5|R8|Aib|Ac|NH2|[A-Z]|[a-z]|[0-9]|' + '|'.join(self.additional_residues.keys()) + ')'
            self.pattern_str = '('+'|'.join(self.additional_residues.keys()) + '|S3|S5|S8|R3|R5|R8|Aib|Ac|NH2|[A-Z]|[a-z]|[0-9])'
        else:
            self.pattern_str = '(S3|S5|S8|R3|R5|R8|Aib|Ac|NH2|[A-Z]|[a-z]|[0-9])'

    def _seq_to_list(self, seq: str) -> list[str]:
        # 需要增加-X-的判断, 比如'RK-M1R-QQ' -> ['R', 'K', 'M1R', 'Q', 'Q']
        # re_aa = re.compile(r'(S3|S5|S8|R3|R5|R8|Aib|Ac|NH2|[A-Z]|[a-z]|[0-9])')
        re_aa = re.compile(r'' + self.pattern_str)
        seq_to_list = re_aa.findall(seq)
        seq_to_list = [x.replace('-', '') for x in seq_to_list if x.strip() != '']
        return [x.replace('X', 'S5') for x in seq_to_list]
    
    def is_stapled(self, seq: str) -> bool:
        '''
            Check if sequence is stapled
        '''
        return len(set(self._seq_to_list(seq)) & set(['X', 'S3', 'S5', 'S8', 'R3', 'R5', 'R8'])) > 0

    def check_seq_validation(self, seq: str) -> None:
        '''
            Check if sequence is valid
        '''
        if self.additional_residues is None:
            additional_residues = {}
        else:
            additional_residues = self.additional_residues

        for step, _s in enumerate(self._seq_to_list(seq)):
            if _s not in self.aa_reversed_dict.keys() and _s not in additional_residues.keys():
                raise ValueError(f'{_s}{step+1} is not a valid amino acid')

    def _one_to_three(self, seq: str) -> str:
        '''
            Convert one letter amino acid to three letter amino acid
        '''
        seq_list = self._seq_to_list(seq)
        # if additional_residues is None:
        self.check_seq_validation(seq)
        
        three_letter_seq = [self.aa_reversed_dict[aa] if aa in self.aa_reversed_dict else aa for aa in seq_list]
        return ' '.join(three_letter_seq)

if __name__ == '__main__':

    seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'
    _task_cid = '88842ae8-d32f-58d7-82d8-046d36f68f42'
    # pathname = os.path.join('data', _task_cid)
    pathname = os.path.join('/home/dahuilangda/Downloads/pep_zju', _task_cid)

    pcp = PhysicochemicalPredictor(sequence=seq, 
                                   topology_file=os.path.join(pathname, 'pep_vac.prmtop'),
                                   trajectory_file=os.path.join(pathname, 'traj.nc'))
    print('helix percent: ', pcp.calc_helix_percent())
    print('sheet percent: ', pcp.calc_extend_percent())
    print('loop percent: ', pcp.calc_loop_percent())
    print('mean bfactor: ', pcp.calc_mean_bfactor())
    print('mol surf: ', pcp.calc_mean_molsurf())
    print('mean gyrate: ', pcp.calc_mean_gyrate())
    print('hydrophobic index: ', pcp.calc_hydrophobic_index(os.path.join(pathname, 'mean_structure.pdb')))
    print('psa: ', pcp.calc_psa(os.path.join(pathname, 'mean_structure.pdb')))
    print('total number of hydrogen bonds: ', pcp.calc_n_hbonds())

    pps = ProtParamsSeq(seq)
    print('length: ', pps.seq_length)
    print('weight: ', pps.weight)
    print('hydrophobicity index: ', pps.hydrophobicity_index)
    print('charge', pps.calc_charge())
    print('charge_density', pps.calc_charge_density())
    print('aromaticity', pps.aromaticity)
    print('fraction_arginine', pps.fraction_arginine)
    print('fraction_lysine', pps.fraction_lysine)
    print('lyticity index: ', pps.lyticity_index)
    print('isoelectric_point: ', pps.isoelectric_point)
