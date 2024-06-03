import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from itertools import combinations
import upsetplot

#-------PREPROCESSING---------------------------------------------------------------------

def extract_first_word_after_dot(filepath):
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    first_word = base_name.split(".")[1] if "." in base_name else base_name
    if first_word =="gprofiler2":
        first_word = "seeds"
    return first_word

def process_file(input_file, output_directory, top_n):
    data = pd.read_csv(input_file, sep="\t")
    top_rows = data.sort_values(by="p_value").head(top_n)
    output_file = os.path.join(output_directory, f"{extract_first_word_after_dot(input_file)}_top_{top_n}_results.tsv")
    top_rows.to_csv(output_file, sep="\t", index=False)
    return top_rows

def adjust_source_names(dataframe):
    dataframe['source'] = dataframe['source'].apply(lambda x: 'GO_ALL' if x.startswith('GO') else x)
    return dataframe

#-------JACCARD-INDEX-----------------------------------------------------------------------------

def calculate_jaccard_index(pathway_sets):
    jaccard_matrix = np.zeros((len(pathway_sets), len(pathway_sets)))
    for i, (file1, pathways1) in enumerate(pathway_sets.items()):
        for j, (file2, pathways2) in enumerate(pathway_sets.items()):
            common_term_ids = pathways1['term_ids'].intersection(pathways2['term_ids'])
            intersection = len(common_term_ids)
            union = len(pathways1['term_ids'].union(pathways2['term_ids']))
            jaccard_index = intersection / union if union != 0 else 0.0
            jaccard_matrix[i, j] = jaccard_index
    return jaccard_matrix

def plot_jaccard_matrix_all(jaccard_matrix, pathway_sets, output_directory, title_prefix=""):
    if np.all(jaccard_matrix == 1) or np.all(jaccard_matrix - np.diag(np.diag(jaccard_matrix)) == 0) or np.all(jaccard_matrix == 0):
        print(f"No data available to plot for {title_prefix}. Because there is no overlap")
        return

    sorted_keys = sorted(pathway_sets.keys())
    sorted_labels = [extract_first_word_after_dot(key) for key in sorted_keys]

    # Reorder the Jaccard matrix to match the sorted labels
    indices = [list(pathway_sets.keys()).index(key) for key in sorted_keys]
    sorted_jaccard_matrix = jaccard_matrix[indices, :][:, indices]

    plt.figure(figsize=(6, 5))
    plt.imshow(sorted_jaccard_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar(label="Jaccard Index")
    plt.clim(0, 1)
    plt.xticks(np.arange(len(sorted_labels)), sorted_labels, rotation=90)
    plt.yticks(np.arange(len(sorted_labels)), sorted_labels)
    plt.title(f"{title_prefix} Jaccard Index Between Top Pathway Sets")
    threshold = (sorted_jaccard_matrix.max() + sorted_jaccard_matrix.min()) / 2
    for i in range(len(sorted_labels)):
        for j in range(len(sorted_labels)):
            color = 'white' if sorted_jaccard_matrix[i, j] < threshold else 'black'
            plt.text(j, i, f"{sorted_jaccard_matrix[i, j]:.2f}", ha='center', va='center', color=color)
    plt.tight_layout()
    plot_file_name = f"{title_prefix.lower().replace(' ', '_')}_jaccard_matrix.png"
    plt.savefig(os.path.join(output_directory, plot_file_name))
    plt.close()

def plot_jaccard_matrix_topX(jaccard_matrix, pathway_sets, output_directory, title_prefix="", top_n=50):
    if np.all(jaccard_matrix == 1) or np.all(jaccard_matrix - np.diag(np.diag(jaccard_matrix)) == 0) or np.all(jaccard_matrix == 0):
        print(f"No data available to plot for {title_prefix}. Because there is no overlap")
        return

    sorted_keys = sorted(pathway_sets.keys())
    sorted_labels = [extract_first_word_after_dot(key) for key in sorted_keys]

    # Reorder the Jaccard matrix to match the sorted labels
    indices = [list(pathway_sets.keys()).index(key) for key in sorted_keys]
    sorted_jaccard_matrix = jaccard_matrix[indices, :][:, indices]

    plt.figure(figsize=(6, 5))
    plt.imshow(sorted_jaccard_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar(label="Jaccard Index")
    plt.clim(0, 1)
    plt.xticks(np.arange(len(sorted_labels)), sorted_labels, rotation=90)
    plt.yticks(np.arange(len(sorted_labels)), sorted_labels)
    plt.title(f"{title_prefix} Jaccard Index Between Top {top_n} Pathway Sets")
    threshold = (sorted_jaccard_matrix.max() + sorted_jaccard_matrix.min()) / 2
    for i in range(len(sorted_labels)):
        for j in range(len(sorted_labels)):
            color = 'white' if sorted_jaccard_matrix[i, j] < threshold else 'black'
            plt.text(j, i, f"{sorted_jaccard_matrix[i, j]:.2f}", ha='center', va='center', color=color)
    plt.tight_layout()
    plot_file_name = f"{title_prefix.lower().replace(' ', '_')}_jaccard_matrix_top_{top_n}.png"
    plt.savefig(os.path.join(output_directory, plot_file_name))
    plt.close()

#-------Overlap-Upset-----------------------------------------------------------------------------

def calculate_overlap_counts(pathway_sets):
    index_tuples = []
    overlap_counts = []

    term_sets = {key: set(data['term_ids']) for key, data in pathway_sets.items()}
    
    all_sets = list(term_sets.keys())
    sorted_keys = sorted(all_sets)

    # Betrachte alle Kombinationen
    for r in range(1, len(term_sets) + 1):
        for subset in combinations(all_sets, r):
            # Berechne die Schnittmenge für diese Kombination
            intersection = set.intersection(*(term_sets[key] for key in subset))
            # Entferne alle Elemente, die in anderen Sets außerhalb dieser Kombination vorkommen
            for other in set(all_sets) - set(subset):
                intersection -= term_sets[other]
            if intersection:  # Ignoriere leere Schnittmengen
                membership = [key in subset for key in sorted_keys]
                index_tuples.append(tuple(membership))
                overlap_counts.append(len(intersection))

    index = pd.MultiIndex.from_tuples(index_tuples, names=sorted_keys)
    overlap_data = pd.Series(overlap_counts, index=index)
    return overlap_data

def plot_overlap_upset(overlap_data, output_directory, title_prefix=""):
    plt.figure(figsize=(15, 20), constrained_layout=True)
    upsetplot.plot(overlap_data, sort_by='cardinality', show_counts=True)
    plt.title(f"UpSet Plot for the Pathway Sets", pad=20, loc ='right')

    plot_file_name = f"{title_prefix.lower().replace(' ', '_')}_overlap_upset_plot.png"
    plt.savefig(f"{output_directory}/{plot_file_name}")
    plt.close()

#-------JACCARD-INDEX-PLOTS-FOR-EACH-SOURCE-----------------------------------------------------------------------------

def create_plots_for_each_source_all(pathway_sets, output_directory):
    unique_sources = set()
    for pathways in pathway_sets.values():
        unique_sources.update(pathways["data"]["source"].unique())
    
    for source in unique_sources:
        source_specific_sets = {}
        for file_key, pathways in pathway_sets.items():
            filtered_data = pathways["data"][pathways["data"]["source"] == source]
            if not filtered_data.empty:
                source_specific_sets[file_key] = {
                    "data": filtered_data,
                    "term_ids": set(filtered_data["term_id"]),
                }
        
        if source_specific_sets:
            jaccard_matrix = calculate_jaccard_index(source_specific_sets)
            plot_jaccard_matrix_all(jaccard_matrix, source_specific_sets, output_directory, f"{source}_")


def create_plots_for_each_source_topX(pathway_sets, output_directory, top_n=50):
    unique_sources = set()
    for pathways in pathway_sets.values():
        unique_sources.update(pathways["data"]["source"].unique())
    
    for source in unique_sources:
        source_specific_sets = {}
        for file_key, pathways in pathway_sets.items():
            filtered_data = pathways["data"][pathways["data"]["source"] == source]
            if not filtered_data.empty:
                source_specific_sets[file_key] = {
                    "data": filtered_data,
                    "term_ids": set(filtered_data["term_id"]),
                }
        
        if source_specific_sets:
            jaccard_matrix = calculate_jaccard_index(source_specific_sets)
            plot_jaccard_matrix_topX(jaccard_matrix, source_specific_sets, output_directory, f"{source}_top_{top_n}", top_n)


#-------VERTEILUNG-DER_QUELLEN-TEXT-UND-PLOT-----------------------------------------------------------------------------


def write_and_plot_source_distribution_normalized(pathway_sets, output_directory):
    output_file = os.path.join(output_directory, "all_pathways_distribution.txt")
    all_source_names = set()
    for data in pathway_sets.values():
        all_source_names.update(data["data"]["source"].unique())
    all_source_names = sorted(list(all_source_names))

    with open(output_file, "w") as f:
        # Header der Datei schreiben
        f.write("File\t" + "\t".join(all_source_names) + "\tTotal Count\n")
        # Bereite Daten für das Plot vor
        plot_data = {}

        for file_key, data in pathway_sets.items():
            total_count = len(data["data"])
            counts = data["data"]["source"].value_counts()
            percentages = [(counts.get(source, 0) / total_count * 100) for source in all_source_names]
            # Berechne prozentuale Anteile für jede Quelle und schreibe sie in die Datei
            row = [file_key] + [f"{p:.2f}%" for p in percentages] + [total_count]
            f.write("\t".join(map(str, row)) + "\n")
            # Speichere Daten für das Plot
            plot_data[file_key] = percentages

    # Erstelle das Balkendiagramm
    fig, ax = plt.subplots(figsize=(6, 2))
    modules = list(plot_data.keys())
    sources = all_source_names
    source_data = np.array([plot_data[module] for module in modules])

    # Erstelle gestapelte Balken für jedes Modul
    bottom = np.zeros(len(modules))
    for i, source in enumerate(sources):
        ax.bar(modules, source_data[:, i], bottom=bottom, label=source)
        bottom += source_data[:, i]

    ax.set_ylabel('Percentage')
    ax.set_title('Percentage Distribution of Sources by Module Normalized')
    ax.legend(title="Sources")

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, 'source_distribution_normalized.png'))
    plt.close()


def write_and_plot_source_distribution(pathway_sets, output_directory):
    output_file = os.path.join(output_directory, "all_pathways_distribution.txt")
    all_source_names = set()
    for data in pathway_sets.values():
        all_source_names.update(data["data"]["source"].unique())
    all_source_names = sorted(list(all_source_names))

    with open(output_file, "w") as f:
        f.write("File\t" + "\t".join(all_source_names) + "\tTotal Count\n")
        plot_data = {}

        for file_key, data in pathway_sets.items():
            total_count = len(data["data"])
            counts = data["data"]["source"].value_counts()
            count_values = [counts.get(source, 0) for source in all_source_names]
            percentages = [(counts.get(source, 0) / total_count * 100) if total_count > 0 else 0 for source in all_source_names]
            row = [file_key] + [f"{p:.2f}%" for p in percentages] + [total_count]
            f.write("\t".join(map(str, row)) + "\n")
            plot_data[file_key] = (count_values, percentages, total_count)

    fig, ax = plt.subplots(figsize=(12, 8))
    modules = list(plot_data.keys())
    bottoms = np.zeros(len(modules))
    last_text_positions = np.zeros(len(modules))

    threshold = 5  # Prozentsatz-Schwelle
    horizontal_offsets = np.linspace(0.1, 0.9, len(all_source_names))  # Generiert relative Positionen innerhalb der Balken

    for i, source in enumerate(all_source_names):
        heights = [plot_data[module][0][i] for module in modules]
        percentages = [plot_data[module][1][i] for module in modules]
        bars = ax.bar(modules, heights, bottom=bottoms, label=source)
        
        for idx, (bar, percentage) in enumerate(zip(bars, percentages)):
            if percentage > 0:  # Nur anzeigen, wenn der Prozentsatz größer als 0 ist
                y_pos = bar.get_y() + bar.get_height() / 2
                x_pos = bar.get_x() + bar.get_width() / 2  # Standardposition in der Mitte des Balkens

                if percentage < threshold:
                    x_pos += bar.get_width() * (horizontal_offsets[i] - 0.5)  # Anpassung nur bei kleinen Werten
                
                ax.text(x_pos, y_pos, f"{percentage:.1f}%", ha='center', va='center')
        
        bottoms += heights

    # Display the total size above each stack, adjusting to avoid overlap with the last percentage
    for j, module in enumerate(modules):
        total_count = plot_data[module][2]
        total_y_pos = bottoms[j]
        if total_y_pos - last_text_positions[j] < 0.1:
            total_y_pos += 0.1
        ax.text(j, total_y_pos, f"Total: {total_count}", ha='center', va='bottom', fontsize=9, color='red')

    ax.set_ylabel('Number of Terms')
    ax.set_title('Distribution and Percentage of Sources by Module')
    ax.legend(title="Sources", loc='upper left', bbox_to_anchor=(1, 1))

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, 'source_distribution_combined.png'))
    plt.close()

#----Textdateien-mit-weiteren-Infos-JI-und-gleiche-Terms-----------------------------

def write_jaccard_index_results_all(pathway_sets, jaccard_matrix, output_directory):
    output_file = os.path.join(output_directory, "jaccard_index_results.tsv")
    labels = [extract_first_word_after_dot(key) for key in pathway_sets.keys()]
    
    with open(output_file, "w") as f:
        f.write("Module 1\tModule 2\tJaccard Index\tCommon Terms Count\n")
        for i, label1 in enumerate(labels):
            for j, label2 in enumerate(labels):
                # Berechnung der gemeinsamen Terme, auch wenn ein Modul mit sich selbst verglichen wird
                common_terms = len(pathway_sets[list(pathway_sets.keys())[i]]['term_ids'].intersection(pathway_sets[list(pathway_sets.keys())[j]]['term_ids']))
                jaccard_index = jaccard_matrix[i, j]
                f.write(f"{label1}\t{label2}\t{jaccard_index}\t{common_terms}\n")

def write_jaccard_index_results_topX(pathway_sets, output_directory):
    output_file = os.path.join(output_directory, "jaccard_with_common_terms.tsv")
    with open(output_file, "w") as f:
        f.write("Module 1\tModule 2\tJaccard Index\tCommon Terms\n")
        for i, (file1, pathways1) in enumerate(pathway_sets.items()):
            for j, (file2, pathways2) in enumerate(pathway_sets.items()):
                common_term_ids = pathways1['term_ids'].intersection(pathways2['term_ids'])
                intersection = len(common_term_ids)
                union = len(pathways1['term_ids'].union(pathways2['term_ids']))
                jaccard_index = intersection / union if union != 0 else 0.0
                common_terms_str = ", ".join(common_term_ids)
                f.write(f"{file1}\t{file2}\t{jaccard_index}\t{common_terms_str}\n")

#-----Aufruf der Methoden-------------------------

def validate(input_directory, output_directory, modus, topX):
    all_pathway_sets = {}
    files = os.listdir(input_directory)
    if modus == "all":
        for file in files:
            if file.endswith(".all_enriched_pathways.tsv"):
                input_file = os.path.join(input_directory, file)
                pathway_data = pd.read_csv(input_file, sep="\t")
                #pathway_data = adjust_source_names(pathway_data) --> war dazu da um alle GO Terms zusammen zu fassen
                file_key = extract_first_word_after_dot(file)
                all_pathway_sets[file_key] = {
                    "data": pathway_data,
                    "term_ids": set(pathway_data["term_id"]),
                }
        #kalkuliere Jaccard Werte und erstelle eine Upset für alle AMIMS mit allen Sources
        jaccard_matrix = calculate_jaccard_index(all_pathway_sets)
        plot_jaccard_matrix_all(jaccard_matrix, all_pathway_sets, output_directory, "Overall ")

        #berechne den totalen overlap und erstelle ein Upset für alle AMIMs mit allen Sources
        overlap_counts = calculate_overlap_counts(all_pathway_sets)
        plot_overlap_upset(overlap_counts, output_directory, "Overall ")

        #kalkuliere Jaccard Werte und erstelle eine Upset für alle AMIMS pro Source
        create_plots_for_each_source_all(all_pathway_sets, output_directory)

        #schreibt die prozentuale Verteilung der Sources in den verschiedenen Amims auf
        #write_source_distribution(all_pathway_sets, output_directory)
        
        write_and_plot_source_distribution(all_pathway_sets, output_directory)
        write_jaccard_index_results_all(all_pathway_sets, jaccard_matrix, output_directory)  
        
    elif modus.startswith("top"):
        for file in files:
            if file.endswith(".all_enriched_pathways.tsv"):
                input_file = os.path.join(input_directory, file)
                pathway_data = process_file(input_file, output_directory, topX)
                #pathway_data = adjust_source_names(pathway_data)
                file_key = extract_first_word_after_dot(file)
                all_pathway_sets[file_key] = {
                    "data": pathway_data,
                    "term_ids": set(pathway_data["term_id"]),
                }
        
        jaccard_matrix = calculate_jaccard_index(all_pathway_sets)
        plot_jaccard_matrix_topX(jaccard_matrix, all_pathway_sets, output_directory, "Overall ", topX)

        overlap_counts = calculate_overlap_counts(all_pathway_sets)
        plot_overlap_upset(overlap_counts, output_directory, "Overall ")
        
        
        create_plots_for_each_source_topX(all_pathway_sets, output_directory, topX)

        #write_source_distribution(all_pathway_sets, output_directory)
    
        write_and_plot_source_distribution(all_pathway_sets, output_directory)
        
        write_jaccard_index_results_topX(all_pathway_sets, output_directory)

    else:
        print("Wrong modus")

#-----Validierung der Module Sizes-----------------------------------------------------------

def tsv_size_plot (tsv_directory, output_directory):
    file_sizes = {}
    for filename in os.listdir(tsv_directory):
        if filename.endswith(".tsv"): 
            key = extract_first_word_after_dot(filename)
            if key != "nodes":  
                with open(os.path.join(tsv_directory, filename), 'r') as file:
                    lines = file.readlines()
                    count = len(lines) - 1  
                    key = extract_first_word_after_dot(filename)
                    file_sizes[key] = count

    if file_sizes:
        labels = list(file_sizes.keys())
        values = list(file_sizes.values())

        plt.figure(figsize=(10, 6))
        bars = plt.bar(labels, values, color='blue')
        plt.xlabel('AMIM')
        plt.ylabel('nodes')
        plt.title('Module Sizes of AMIMs')
        plt.xticks(rotation=45)
        plt.tight_layout() 
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), ha='center', va='bottom')

        plt.show()

    plt.savefig(os.path.join(output_directory, 'module_sizes.png'))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Execute validation")

    parser.add_argument("--input_dir", required=True, help="Input directory path")
    parser.add_argument("--output_dir", required=True, help="output directory path")
    parser.add_argument("--modus", required=True, help="validation for all or topX")
    parser.add_argument("--topX", required=False, help=" topX must be an Integer, for number of terms per file that should be invested")
    parser.add_argument("--tsv_dir", required=True, help="tsv directory of the modules generated from the AMIMs")
    
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    modus = args.modus
    topX = args.topX
    if topX is not None:
        topX = int(topX)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    validate(input_dir, output_dir, modus, topX)
    tsv_size_plot(args.tsv_dir, output_dir)