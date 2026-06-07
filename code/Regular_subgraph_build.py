import networkx as nx
import statsmodels.api as sm
import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
#############
TF_CANDIDATE_TOP_N = 20


def analyze_gexf_graph(gexf_file_path, output_file):
    """
    Analyzes a GEXF graph and saves cycle information to a TSV file.

    Args:
        gexf_file_path (str): Path to the GEXF file.
        output_file (str): Path to save the TSV file.

    Returns:
        None
    """

    def read_gexf_file(file_path):
        return nx.read_gexf(file_path)

    def find_cycles(graph):
        return nx.simple_cycles(graph)

    def get_cycle_info(cycles):
        cycle_info = []
        for cycle in cycles:
            num_nodes = len(cycle)
            num_edges = sum(graph[cycle[i]][cycle[i + 1]]['weight'] for i in range(num_nodes - 1))
            cycle_info.append((num_nodes, num_edges))
        return cycle_info

    def save_to_tsv(cycle_info, output_file):
        df = pd.DataFrame(cycle_info, columns=['Nodes', 'Edges'])
        df.to_csv(output_file, sep='\t', index=False)

    graph = read_gexf_file(gexf_file_path)
    cycles = find_cycles(graph)
    cycle_info = get_cycle_info(cycles)
    save_to_tsv(cycle_info, output_file)

    print(f" {len(cycles)} ， {output_file} 。")


def read_gexf_file(file_path):
    return nx.read_gexf(file_path)
def is_weighted_graph(graph):
    for u, v, data in graph.edges(data=True):
        if 'weight' in data:
            return True  #  'weight' ，
    return False  # ，

def build_subgraph(graph, target_node):
    visited = set()
    subgraph = nx.DiGraph()

    queue = [target_node]
    visited.add(target_node)

    while queue:
        current_node = queue.pop(0)
        for neighbor in graph.predecessors(current_node):
            if neighbor not in visited:
                queue.append(neighbor)
                visited.add(neighbor)
                subgraph.add_edge(neighbor, current_node, weight=graph[neighbor][current_node]['weight'])

    return subgraph


def build_subgraph_cycle(graph, target_node):
    visited = set()
    subgraph = nx.DiGraph()

    queue = [target_node]
    visited.add(target_node)

    while queue:
        current_node = queue.pop(0)
        for neighbor in graph.predecessors(current_node):
            if neighbor not in visited:
                queue.append(neighbor)
                visited.add(neighbor)
                subgraph.add_edge(neighbor, current_node, weight=graph[neighbor][current_node]['weight'])

    for node in subgraph.nodes():
        for neighbor in graph.predecessors(node):
            if neighbor in subgraph.nodes():
                subgraph.add_edge(neighbor, node, weight=graph[neighbor][node]['weight'])

    return subgraph

def build_subgraph_cycle_no_weight(graph, target_node):
    visited = set()
    subgraph = nx.DiGraph()

    queue = [target_node]
    visited.add(target_node)

    while queue:
        current_node = queue.pop(0)
        for neighbor in graph.predecessors(current_node):
            if neighbor not in visited:
                queue.append(neighbor)
                visited.add(neighbor)
                subgraph.add_edge(neighbor, current_node)

    for node in subgraph.nodes():
        for neighbor in graph.predecessors(node):
            if neighbor in subgraph.nodes() and not subgraph.has_edge(neighbor, node):
                subgraph.add_edge(neighbor, node)

    return subgraph

def build_subgraph_from_gexf_plt(graph, target_node):
    target_node1 = target_node

    if target_node1 not in graph:
        print(f"The node {target_node1} is not in the digraph.")
        return None  # : raise nx.NetworkXError(f"The node {target_node} is not in the digraph.")
    subgraph = build_subgraph_cycle(graph, target_node1)

    #for edge in subgraph.edges(data=True):
        #print(f"Edge from {edge[0]} to {edge[1]}, weight: {edge[2]['weight']}")

    num_nodes = subgraph.number_of_nodes()
    num_edges = subgraph.number_of_edges()

    print(f": {num_nodes}")
    print(f": {num_edges}")

    return subgraph

def build_subgraph_from_gexf_noweight(graph, target_node):
    target_node1 = target_node
    if target_node1 not in graph:
        print(f"The node {target_node1} is not in the digraph.")
        return None
    subgraph = build_subgraph_cycle_no_weight(graph, target_node1)
    if not subgraph.nodes():
        print(f"No predecessors found for {target_node1}, returning empty subgraph.")
        return None

    num_nodes = subgraph.number_of_nodes()
    num_edges = subgraph.number_of_edges()
    print(f"Subgraph nodes: {num_nodes}, edges: {num_edges}")
    return subgraph

def build_subgraph_from_gexf(gexf_path, target_node):
    graph = read_gexf_file(gexf_path)
    subgraph = build_subgraph(graph, target_node)

    num_nodes = subgraph.number_of_nodes()
    num_edges = subgraph.number_of_edges()

    print(f": {num_nodes}")
    print(f": {num_edges}")

    return subgraph


def normalize_graph(graph0):
    for node in graph0.nodes():
        total_weight = sum(graph0[node][neighbor]['weight'] for neighbor in graph0.successors(node))
        for neighbor in graph0.successors(node):
            if total_weight > 0:
                graph0[node][neighbor]['weight'] = graph0[node][neighbor]['weight'] / total_weight
#######################################################################################
#######################################################################################
def katz_influence_on_target(graph, target, alpha=0.1, beta=1.0):
    target = target[0]
    katz_centrality = nx.katz_centrality(graph, alpha=alpha, beta=beta, weight='weight')

    influence = {node: 0 for node in graph.nodes() if node != target}

    for node in graph.nodes():
        if node != target:
            paths = nx.all_simple_paths(graph, source=node, target=target)

            for path in paths:
                path_length = len(path) - 1  # 
                influence[node] += (alpha ** path_length) * katz_centrality[node]



    return influence

def betweenness_centrality_influence(graph, target):
    target = target[0]
    betweenness_centrality = nx.betweenness_centrality(graph, weight='weight', normalized=True)

    influence = {node: betweenness_centrality[node] for node in graph.nodes() if node != target}

    #sorted_influence = sorted(influence.items(), key=lambda item: item[1], reverse=True)

    return influence

def betweenness_centrality_influence4target(graph, target):
    target = target[0]
    influence = {node: 0 for node in graph.nodes() if node != target}

    for source in graph.nodes():
        if source == target:
            continue  # 

        try:
            paths = list(nx.all_shortest_paths(graph, source=source, target=target, weight='weight'))
        except nx.NetworkXNoPath:
            continue  # ，

        for path in paths:
            for node in path[1:-1]:
                influence[node] += 1

    return influence


def linear_threshold_influence(graph, target):
    influence = {node: 0 for node in graph.nodes() if node != target}
    direct_neighbors = [] #

    for node in graph.nodes():
        if graph.has_edge(node, target):
            direct_neighbors.append(node)
            if 'weight' in graph[node][target]:
                influence[node] = graph[node][target]['weight']
            else:
                influence[node] = 1  # ，1

    for node in graph.nodes():
        if node != target:
            for neighbor in graph.successors(node):
                if graph.has_edge(node, neighbor) and neighbor != target:
                    if 'weight' in graph[node][neighbor]:
                        influence[node] += (graph[node][neighbor]['weight'] * influence[neighbor])
                    else:
                        influence[node] += influence[neighbor] #，weight1

    for node in direct_neighbors:
        influence[node] += 1000000

    return influence

def restarted_random_walk_influence(G, target_node, num_walks=10000, walk_length=5):
    target_node = target_node[0]
    if not G.nodes():  # 
        print(f"Subgraph for {target_node} is empty, skipping random walk.")
        return {node: 0 for node in G.nodes()}  # 
    influence_scores = {node: 0 for node in G.nodes()}

    for walk in range(num_walks):
        current_node = np.random.choice(list(G.nodes()))
        path = [current_node]

        for step in range(walk_length):
            if current_node == target_node:
                influence_scores[path[0]] += 1
                break
            successors = list(G.successors(current_node))
            if successors:
                current_node = np.random.choice(successors)  # 
                path.append(current_node)
            else:
                break

        if current_node != target_node:
            continue

        for node in path:
            influence_scores[node] += 1

    return influence_scores

def weighted_restarted_random_walk_influence(G, target_node, num_walks=10000, walk_length=5):

    target_node = target_node[0]
    influence_scores = {node: 0 for node in G.nodes()}

    for walk in range(num_walks):
        current_node = np.random.choice(list(G.nodes()))
        path = [current_node]  # 

        for step in range(walk_length):
            if current_node == target_node:
                influence_scores[path[0]] += 1
                break
            successors = list(G.successors(current_node))
            if successors:
                weights = [G[current_node][succ]['weight'] for succ in successors]
                total_weight = sum(weights)
                if total_weight > 0:  # 
                    probabilities = [w / total_weight for w in weights]
                    current_node = np.random.choice(successors, p=probabilities)
                    path.append(current_node)  # 
                else:
                    break
            else:
                break

        if current_node != target_node:
            continue

        for node in path:
            influence_scores[node] += 1

    return influence_scores

def get_directly_connected_nodes(graph, target_node):
    #target_node = target_node[0]
    predecessors = list(graph.predecessors(target_node))


    return predecessors

def Independent_Cascade_Model_influence(G, target_node, num_walks=100000, walk_length=20):
    influence_scores = {node: 0 for node in G.nodes()}

    for _ in range(num_walks):
        current_node = np.random.choice(list(G.nodes()))
        for _ in range(walk_length):
            if current_node == target_node:
                influence_scores[current_node] += 1
                break
            successors = list(G.successors(current_node))
            if successors:
                weights = [G[current_node][succ].get('weight', 1) for succ in successors]
                if weights:
                    probabilities = weights / np.sum(weights)
                    current_node = np.random.choice(successors, p=probabilities)
                else:
                    break
            else:
                break

    return influence_scores

def calculate_influence_by_paths(graph, target):
    target = target[0]
    influence = {node: 0 for node in graph.nodes() if node != target}

    for node in graph.nodes():
        if node == target:
            continue
        try:
            shortest_paths = nx.all_shortest_paths(graph, source=node, target=target)
            influence[node] = len(list(shortest_paths))  # 
        except nx.NetworkXNoPath:
            influence[node] = 0  # ， 0

    return influence

#######################################################################################
#######################################################################################


def gene_name_2_ID(bed_file_path, gene_names):
    df = pd.read_csv(bed_file_path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type'])

    name_2_id_map = dict(zip(df['gene_name'], df['gene_id']))

    valid_ids = []
    valid_names = []  # 

    for name in gene_names:
        if name in name_2_id_map:
            valid_ids.append(name_2_id_map[name])
            valid_names.append(name)

    return np.array(valid_ids), np.array(valid_names)



def get_influence_tfgeneID_LTM(graph, target_node=None, numbers=None, expression_TF_gene=None, bed_path=None):


    if is_weighted_graph(graph):
        print("The graph is weighted. Using build_subgraph_from_gexf_plt.")
        sub_graph = build_subgraph_from_gexf_plt(graph, target_node)
    else:
        print("The graph is unweighted. Using build_subgraph_from_gexf_noweight.")
        sub_graph = build_subgraph_from_gexf_noweight(graph, target_node)

    if sub_graph is None:
        return None
    pankrang = linear_threshold_influence(sub_graph, target_node)
    sorted_scores = sorted(pankrang.items(), key=lambda item: item[1], reverse=True)
    top_influencers = sorted_scores[:numbers]
    gene_names = [gene[0] for gene in top_influencers]
    geneNAME = np.array(gene_names, dtype=object).ravel()
    geneNAME = np.intersect1d(geneNAME, expression_TF_gene)
    if bed_path is None:
        bed_path = '/data/lab/wangshixian/GRNTWAS_STAR/data/anno_info/gene.bed'
    tf_gene_ID, _ = gene_name_2_ID(bed_path, geneNAME)
    return tf_gene_ID


def get_influence_tfgeneID_KATZ(gexf_path='/data/lab/wangshixian/GRNTWAS_STAR/data/GRN_data/brain_tf_network.gexf',
                               target_node=None, numbers=None, expression_TF_gene=None):

    sub_graph = build_subgraph_from_gexf_plt(gexf_path, target_node)
    if sub_graph is None:
        return None  # ， None

    influence_scores = katz_influence_on_target(sub_graph, target_node)

    if influence_scores: # influence_scores
        sorted_influencers = sorted(influence_scores.items(), key=lambda item: item[1], reverse=True)
        top_influencers = [item[0] for item in sorted_influencers[:numbers]] # 
    else:
        return np.array([]) # influence_scores，

    gene_names = [gene[0] for gene in top_influencers]
    geneNAME = np.array(gene_names, dtype=object).ravel()
    geneNAME = np.intersect1d(geneNAME, expression_TF_gene)
    tf_gene_ids = gene_name_2_ID('/data/lab/wangshixian/GRNTWAS_STAR/data/anno_info/gene.bed', geneNAME)

    return tf_gene_ids


def get_influence_tfgeneID_paths(sub_graph, target_node=None, numbers=None):
    pankrang = calculate_influence_by_paths(sub_graph, target_node)
    sorted_scores = sorted(pankrang.items(), key=lambda item: item[1], reverse=True)
    top_influencers = sorted_scores[:numbers]
    gene_names = [gene[0] for gene in top_influencers]
    geneNAME = np.array(gene_names, dtype=object).ravel()
    tf_gene_ID = gene_name_2_ID('/data/lab/wangshixian/GRNTWAS_STAR/data/anno_info/gene.bed', geneNAME)
    return tf_gene_ID

def get_influence_tfgeneID_predecessors(sub_graph, target_node=None, numbers=None, expression_TF_gene=None):
    pankrang = get_directly_connected_nodes(sub_graph, target_node)
    gene_names = [gene for gene in pankrang]
    geneNAME = np.array(gene_names, dtype=object).ravel()
    geneNAME = np.intersect1d(geneNAME, expression_TF_gene)
    tf_gene_ID = gene_name_2_ID('/data/lab/wangshixian/GRNTWAS_STAR/data/anno_info/gene.bed', geneNAME)
    return tf_gene_ID

def get_influence_tfgeneID_RRW(gexf_path = '/data/lab/wangshixian/GRNTWAS_STAR/data/GRN_data/brain_tf_network.gexf', target_node=None,numbers=None, expression_TF_gene=None):
    sub_graph = build_subgraph_from_gexf_noweight(gexf_path, target_node)
    if sub_graph is None:
        print(f"No valid subgraph for {target_node}, skipping.")
        return None
    pankrang = restarted_random_walk_influence(sub_graph, target_node)
    sorted_scores = sorted(pankrang.items(), key=lambda item: item[1], reverse=True)
    top_influencers = sorted_scores[:numbers]
    gene_names = [gene[0] for gene in top_influencers]
    print(gene_names)
    geneNAME = np.array(gene_names, dtype=object).ravel()
    geneNAME = np.intersect1d(geneNAME, expression_TF_gene)
    tf_gene_ID = gene_name_2_ID('/data/lab/wangshixian/GRNTWAS_STAR/data/anno_info/gene.bed', geneNAME)
    return tf_gene_ID


def select_TFs_via_graph_guided_bayesian(gexf_path, target_node_name, target_expr_values,
                                         gene_exp_df, sampleID, bed_path, expression_TF_gene_names):

    # -------------------------------------------------------
    # -------------------------------------------------------
    sub_graph = build_subgraph_from_gexf_noweight(gexf_path, target_node_name)
    if sub_graph is None or len(sub_graph.nodes) == 0:
        return []

    try:
        rwr_scores_dict = nx.pagerank(sub_graph.reverse(), alpha=0.85, personalization={target_node_name: 1})
    except Exception as e:
        print(f"RWR Error: {e}")
        return []

    # -------------------------------------------------------
    # -------------------------------------------------------
    candidate_names = [
        name for name in rwr_scores_dict.keys()
        if name in expression_TF_gene_names and name != target_node_name
    ]

    if len(candidate_names) < 2:
        return []

    candidate_names = sorted(candidate_names, key=lambda x: rwr_scores_dict[x], reverse=True)[:TF_CANDIDATE_TOP_N]

    tf_ids, valid_tf_names = gene_name_2_ID(bed_path, candidate_names)

    if len(tf_ids) == 0:
        return []

    final_scores = np.array([rwr_scores_dict[name] for name in valid_tf_names])

    # -------------------------------------------------------
    # -------------------------------------------------------
    y = target_expr_values - np.mean(target_expr_values)

    tf_data = gene_exp_df[gene_exp_df['TargetID'].isin(tf_ids)].copy()
    tf_data = tf_data.drop_duplicates(subset=['TargetID'])  #  ID
    tf_data = tf_data.set_index('TargetID')
    tf_data = tf_data.reindex(tf_ids)  # 

    if tf_data.isnull().values.any():
        valid_mask = ~tf_data.isnull().any(axis=1)
        tf_data = tf_data[valid_mask]
        tf_ids = tf_ids[valid_mask]  #  ID 
        final_scores = final_scores[valid_mask]  # 

    if tf_data.empty:
        return []

    X = tf_data[sampleID].values.T

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    # -------------------------------------------------------
    # -------------------------------------------------------
    epsilon = 1e-4
    max_s = np.max(final_scores)
    normalized_scores = final_scores / (max_s + epsilon)

    penalty_weights = 1.0 / (normalized_scores + 0.05)  # 0.05 ，
    # -------------------------------------------------------
    # -------------------------------------------------------
    X_weighted = X / penalty_weights[np.newaxis, :]

    try:
        model = LassoCV(cv=5, n_jobs=1, random_state=42, max_iter=2000).fit(X_weighted, y)
        beta = model.coef_
        support_indices = np.where(model.coef_ != 0)[0]

        if len(support_indices) == 0:
            return []


    except Exception as e:
        print(f"Lasso step failed: {e}")
        return []

    # -------------------------------------------------------
    # -------------------------------------------------------
    X_subset = X[:, support_indices]
    selected_ids_subset = tf_ids[support_indices]

    try:
        X_subset_with_const = sm.add_constant(X_subset)
        ols_model = sm.OLS(y, X_subset_with_const).fit()
        r2 = ols_model.rsquared_adj
        if ols_model.rsquared_adj < 0.05:
            # print(f"  Target {target_node_name}: Low R2 ({ols_model.rsquared_adj:.4f}), skipped.")
            return []
        p_values = ols_model.pvalues[1:]
        significant_mask = p_values < 0.05

        final_selected_tf_ids = selected_ids_subset[significant_mask].tolist()

        # print(f"  Step 1 Lasso: {len(support_indices)} -> Step 2 OLS: {len(final_selected_tf_ids)}")
        return final_selected_tf_ids

    except Exception as e:
        print(f"OLS step failed: {e}")
        return tf_ids[support_indices].tolist()
