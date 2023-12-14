from glob import glob
from math import ceil

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns


def round_val(val, digit=3):
    return ceil(val * 10 ** digit) / 10 ** digit


def draw_emission(model, save_name):
    sns.set(font_scale=0.5)
    # select rows contains emissionprobs и active == 1.0
    emission = model[(model['type'] == "emissionprobs") & (model['active'] == 1.0)]
    # Getting label names
    marks = emission['mark'].unique()
    # Fill dictionary
    emissions = {"State (Emission order)": [i + 1 for i in range(int(emission['state_from'].unique().max()))]}
    for i in marks:
        emissions.update({i: list(emission[emission['mark'] == i].sort_values(by=['state_from'])['score'])})
    # Create Dataframe from dictionary with index column
    table = pd.DataFrame(emissions).set_index("State (Emission order)")
    # round values to 3 digits
    table = table.applymap(round_val)
    # Draw heatmap without bar
    sns.heatmap(table, cmap="Blues", annot=True, cbar=False, square=True, annot_kws={"fontsize": 3}, linewidths=1)
    # Save heatmap
    plt.savefig(f"{save_name}_emission.png", dpi=800, bbox_inches='tight')
    plt.close()


def draw_transition(model, save_name, show_transition=False):
    size = int(model.iloc[0]['type'])
    # select rows contains transitionprobs
    transition = model[(model['type'] == "transitionprobs")]
    array = [[round_val(float(
        transition[(transition['state_from'] == i + 1) & (transition['state_to'] == str(j + 1))].iloc[0]['mark'])) for
        j in range(size)] for i in range(size)]

    ticks = [i for i in range(1, size + 1)]
    sns.heatmap(
        array,
        xticklabels=ticks,
        yticklabels=ticks,
        cmap="Blues",
        annot=True,
        cbar=False,
        square=True,
        annot_kws={"fontsize": 3},
        linewidths=1
    )
    plt.savefig(f"{save_name}_transition.png", dpi=800, bbox_inches='tight')
    plt.close()


pathes = glob(r"ELIM_MODEL_PATH\elim*.txt")
for path in pathes:
    model = pd.read_csv(path, sep='\t',
                        names={
                            'type': str,
                            'state_from': int,
                            'state_to': int,
                            'mark': str,
                            'active': int,
                            'score': float},
                        )
    draw_transition(model, path.replace(".txt", ''), show_transition=True)
    draw_emission(model, path.replace(".txt", ''))
