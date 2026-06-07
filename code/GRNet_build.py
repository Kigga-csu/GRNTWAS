import pandas as pd
import argparse
import networkx as nx




def extract_filename_from_path(address):
    filename_with_extension = address.split('/')[-1]
    
    filename_without_extension = filename_with_extension.split('.')[0]
    
    return filename_without_extension

def csv_information(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')  # sep
    except Exception as e:
        print(f"An error occurred: {e}")
        exit()

    print(f": {df.shape[0]}")
    print(f": {df.shape[1]}")

    for column in df.columns:
        unique_count = df[column].nunique()
        print(f" '{column}'  {unique_count} ")

    print(":", df.columns.tolist())
    return df

def GRN_extract_1(df, gene_name):
    filtered_data = df[df['gene'] == gene_name]
    print(filtered_data[['TF', 'bestMotif', 'NES', 'Genie3Weight', 'Confidence']])
    return filtered_data

def GRN_build_no_relationship(df):
    df = df.drop_duplicates()
    G = nx.DiGraph()
    for index, row in df.iterrows():
        source = row[0]  # 
        target = row[1]  # 

        G.add_edge(source, target)

    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    print(f": {num_nodes}")
    print(f": {num_edges}")
    return G

def GRN_DLG_build(df, filename):
    G = nx.DiGraph()

    for index, row in df.iterrows():
        G.add_edge(row['TF'], row['gene'], weight=row['NES'])

    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    print(f": {num_nodes}")
    print(f": {num_edges}")

    nx.write_gexf(G, '../data/GRN_data/'+filename+'.gexf')
    return G

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TSV")
    parser.add_argument('file_path', type=str, help="TSV")

    args = parser.parse_args()
    df = csv_information(args.file_path)  # 
    filename = extract_filename_from_path(args.file_path)
    G = GRN_DLG_build(df, filename)  # DataFrame
