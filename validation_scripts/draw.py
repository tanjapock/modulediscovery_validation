import graph_tool.all as gt
import csv
import pandas as pd

def read_seeds(seed_file):
    seeds = set()
    with open(seed_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:  # Sicherstellen, dass die Zeile nicht leer ist
                seeds.add(row[0])
    return seeds

def visualize_graph(file_path, seed_file):
    # Lade den Graphen
    g = gt.load_graph(file_path)

   
    genes_df = pd.read_csv(seed_file, sep='\t')
    if 'gene_id' in genes_df.columns:
        seed_ids = set(genes_df['gene_id'])
    else:
        genes_df = pd.read_csv(seed_file, sep='\t', header=None)
        seed_ids = set(genes_df.iloc[:, 0]) 

    # Erstelle eine neue Vertex-Eigenschaft 'color'
    vertex_color = g.new_vertex_property("vector<double>")
    vertex_text = g.new_vertex_property("string")
    

    # Gehe durch alle Knoten und färbe die Seed-Gene
    for v in g.vertices():
        gene_id = g.vertex_properties["name"][v]
        if g.vertex_properties["name"][v] in seed_ids:
            vertex_color[v] = [1.0, 0.0, 0.0, 1.0]  # Rot für Seed-Gene
            vertex_text[v] = ""
        else:
            vertex_color[v] = [0.0, 0.0, 1.0, 1.0]  # Blau für andere Gene
            vertex_text[v] = gene_id

    # Position der Knoten mit sfdp-Layout
    pos = gt.sfdp_layout(g)

    # Zeichne den Graphen mit den eingefärbten Knoten
    gt.graph_draw(g, pos, vertex_fill_color=vertex_color, output_size=(1000, 1000), vertex_text=vertex_text, vertex_font_size = 12, edge_pen_width=3, output = "thyroid_cancer_16.05_graphs/thyroid_cancer_16.05_robust_d24.png")
    #gt.graph_draw(g, pos, vertex_fill_color=vertex_color, output_size=(1000, 1000),  vertex_font_size = 12, edge_pen_width=1.5, output = "thyroid_cancer_10.05_graphs_new/thyroid_cancer_10.05_robust_d24.png")


if __name__ == "__main__":
    file_path = "/nfs/home/students/t.pock/modulediscovery_2.0/thyroid_cancer_16.05_diamond24/moduleparser/cancer_uniprot_seeds_disgenet.robust.gt"
    seed_file = "/nfs/home/students/t.pock/modulediscovery_2.0/input_files/cancer_uniprot_seeds_disgenet.csv"
    visualize_graph(file_path, seed_file)