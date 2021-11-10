import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp

def load_network(file_path):
    ppi = pd.read_table(filepath_or_buffer=file_path, header=None,
                        index_col=None, names=['source', 'target'], sep='\t')
    ppi_nodes = pd.concat([ppi['source'], ppi['target']], ignore_index=True)
    ppi_nodes = pd.DataFrame(ppi_nodes, columns=['nodes']).drop_duplicates()
    ppi_nodes.reset_index(drop=True, inplace=True)
    return ppi,ppi_nodes

def intersec_net_x_y(netx, nety, nodesx, nodesy):
    nodes_xy = pd.merge(left=nodesx, right=nodesy, on='nodes', how='inner')
    merg1 = pd.merge(left=netx, right=nodes_xy, left_on='source', right_on='nodes', how='left')
    merg1.dropna(axis=0, how='any', inplace=True)
    merg2 = pd.merge(left=merg1, right=nodes_xy, left_on='target', right_on='nodes', how='left')
    merg2.dropna(axis=0, how='any', inplace=True)
    merg2.reset_index(drop=True, inplace=True)
    netx_xy = merg2.loc[:, ['source', 'target']]

    merg1 = pd.merge(left=nety, right=nodes_xy, left_on='source', right_on='nodes', how='left')
    merg1.dropna(axis=0, how='any', inplace=True)
    merg2 = pd.merge(left=merg1, right=nodes_xy, left_on='target', right_on='nodes', how='left')
    merg2.dropna(axis=0, how='any', inplace=True)
    merg2.reset_index(drop=True, inplace=True)
    nety_xy = merg2.loc[:, ['source', 'target']]

    return nodes_xy, netx_xy, nety_xy

def convert_adj_to_edgeset(adj_s_k, node_df):
    sour_lst = []
    targ_lst = []
    for i in np.arange(0, node_df.shape[0]):
        row = adj_s_k[i, i:]
        row = np.pad(row, (i, 0), 'constant', constant_values=(0, 0))
        targ = node_df.loc[row.nonzero()[0], :]['nodes'].tolist()
        sour = [node_df.iloc[i, 0] for j in range(0, len(targ))]
        targ_lst = targ_lst + targ
        sour_lst = sour_lst + sour
    edge_df = pd.DataFrame({'source': sour_lst, 'target': targ_lst})

    return edge_df

def get_nodeset_from_edgeset(edgeset):
    net_1 = edgeset.drop_duplicates('source', 'first', inplace=False)
    net_2 = edgeset.drop_duplicates('target', 'first', inplace=False)
    nodes_df = pd.concat([net_1['source'], net_2['target']],
                        ignore_index=True)
    nodes_df = pd.DataFrame(nodes_df, columns=['nodes']).drop_duplicates('nodes')
    nodes_df.reset_index(drop=True,inplace=True)
    return nodes_df

def calc_gene_weight(mut_df, gene_df ):
    M_lst = []
    mut_num = mut_df.sum(axis=0)
    N_max_recip = 1 / (mut_num.max() + np.spacing(1))

    for i in np.arange(0, gene_df.shape[0]):
        g = gene_df.iloc[i, 0]
        Mg = 0
        if g in mut_df.index:
            Ki = mut_df.loc[g, :]
            Ki = Ki[Ki.to_numpy().nonzero()[0]]
            if len(Ki) > 0:
                for k in Ki.index.tolist():
                    Mg = Mg + 1 / mut_num[k]
            else:
                Mg = N_max_recip
        else:
            Mg = N_max_recip
        M_lst.append(Mg)

    M_recip_lst = [1 / i for i in M_lst]
    M_recip_norm_lst=[]
    arr = np.asarray(M_recip_lst)
    for x in arr:
        M_recip_norm_lst.append(float(x - np.min(arr)) / (np.max(arr) - np.min(arr)))
    return pd.DataFrame({'nodes':gene_df.iloc[:,0].values.tolist(),'weight':M_recip_norm_lst})


def find_penalty_weighted_MCV_single_layer(sample_net, weig_nodes,lamd):
    N1 = weig_nodes.shape[0]
    N2 = sample_net.shape[0]
    A = np.zeros([N2, N1])

    for i in np.arange(0, N2):
        s = sample_net.loc[i, 'source']
        t = sample_net.loc[i, 'target']
        A[i, weig_nodes[weig_nodes.nodes == s].index[0]] = 1
        A[i, weig_nodes[weig_nodes.nodes == t].index[0]] = 1

    m = gp.Model('matrixs1')
    x = m.addMVar(shape=N1, vtype=GRB.BINARY, name='x')
    obj1 = np.ones([1, N1])
    obj = weig_nodes['weight'].to_numpy().reshape([1,N1])
    m.setObjective(obj1 @ x + lamd * obj @ x , GRB.MINIMIZE)
    A = sp.csr_matrix(A)
    rhs = np.ones([1, N2])[0]
    m.addConstr(A @ x >= rhs, name='c')
    m.setParam('OutputFlag', 0)
    m.optimize()
    drivers = weig_nodes[x.X == 1]
    return drivers['nodes'].values.tolist()
