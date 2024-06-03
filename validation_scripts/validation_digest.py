#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_p_values(folder_path):
    # Prüfen, ob der Pfad existiert
    if not os.path.exists(folder_path):
        print("Der angegebene Pfad existiert nicht.")
        return

    # Vorbereitung der Plotdaten
    plot_data = {}
    #colors = plt.cm.get_cmap('Accent', 4)  # Farbkarte für bis zu 10 verschiedene Dateien
    color_map = {
        
        'diamond': 'green',
        'domino': 'orange',
        'seeds': 'purple',
        'robust': 'blue'
    }

    # Durchsuchen aller Unterordner und Dateien
    for subdir, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('p-value_validation.csv'):
                file_path = os.path.join(subdir, file)
                # Bestimmung des Namens für die Legende
                subdir_name = os.path.basename(subdir)
                if '.' in subdir_name:
                    plot_name = subdir_name.split('.')[1]
                else:
                    plot_name = 'seeds'
                
                # Daten einlesen
                data = pd.read_csv(file_path, skiprows=0)
                data.columns = ['Source', 'pvalue']
                
                # -log10 der p-Werte berechnen
                data['-log10(pvalue)'] = -np.log10(data['pvalue'])
                
                # Daten zum Plotten hinzufügen
                plot_data[plot_name] = data

    sorted_plot_data = dict(sorted(plot_data.items()))
      
    # Erstellen des Plots
    plt.figure(figsize=(6, 4))
    for key, value in sorted_plot_data.items():
        plt.scatter(value['Source'], value['-log10(pvalue)'], color=color_map.get(key, 'black'), label=key)

    # Threshold-Linie einzeichnen
    plt.axhline(-np.log10(0.05), color='red', linestyle='dashed', linewidth=1, label='Threshold')

    plt.xlabel('Source')
    plt.ylabel('-log10(pvalue)')
    plt.title('P-Values for Functional Coherence ')
    handles, labels = plt.gca().get_legend_handles_labels()
    sorted_labels_handles = sorted(zip(labels, handles), key=lambda x: x[0])
    labels, handles = zip(*sorted_labels_handles)

    plt.legend(handles, labels)
    plt.legend()
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 3.5)
    plt.tight_layout()


    # Speichern des Plots als PNG-Datei
    plt.savefig("/nfs/home/students/t.pock/modulediscovery_2.0/thyroid_cancer_10.05_validation_digest/pvalue_validation_new.png", format='png')

    # Anzeigen des Plots
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Bitte geben Sie den Pfad zum Ordner als Argument an.")
    else:
        plot_p_values(sys.argv[1])
