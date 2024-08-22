import io
import os
import sys
import time
import math
import argparse
import requests
import json
import vcf

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
from adjustText import adjust_text

# for seaborn plots, change size of canvas and create color palette object
sns.set_theme(rc={'figure.figsize': (11.7, 8.27)})
color_palette = None


# required installations
# pip install graphviz
# pip install pydot
os.environ["PATH"] += r"C:\goutam\soham\Python\Graphviz\bin" # local path to Graphviz/bin

# initialize user input parser and arguments
parser = argparse.ArgumentParser(description='Foo')
parser.add_argument('--vcf_dir', metavar='vcf_dir', type=str, required=True)
parser.add_argument('--vcf_suff', metavar='vcf_suff', type=str, required=True)
parser.add_argument('--tree_dir', metavar='tree_dir', type=str, required=True)
parser.add_argument('--neoag_dir', metavar='neoag_dir', type=str, required=True)
parser.add_argument('--patient', metavar='patient', type=str, required=False)
args = parser.parse_args()

# convert arguments into global variables
vcf_directory = args.vcf_dir
patient = args.patient
vcf_suffix = args.vcf_suff
neoant_directory = args.neoag_dir
tree_directory = args.tree_dir

# create patient, mutation, and clone summary dataframes
patient_summary = pd.DataFrame(columns=['Patient', 'Sample'], index=[])
mutation_summary = pd.DataFrame()
clone_summary = pd.DataFrame()

# create vcf, neoantigen, and tree dictionaries to organize individual datasets
vcf_dict = {}
neoantigen_dict = {}
tree_dict = {}

# graph_df object for storing temporary dataframes used in plotting graphs and drawing trees
graph_df = None

# output folder
output_path = "C:/goutam/soham/MSKCC/summer_2024/"


# insert value(s) into dataframe column(s)
def insert_value(df, index, cols, vals):
    for i in range(0, len(cols)):
        df.loc[index, cols[i]] = vals[i]


# return subset of df where column matches value
def filter_df_by_value(df, col, val, check_substring):
    if check_substring:
        return df[df[col].str.contains(val)]
    else:
        return df[df[col] == val]


# replace occurrence of value with new value in dataframe column
def replace_in_col(col, val, new_val):
    col.replace(val, new_val, inplace=True)


# splice given string based on specific start and end points
def splice(string):
    start = string.find("|") + 1
    end = string[start:].find("|") + start
    return string[start:end]


# create pdf object, and the axes that will need to plotted and adjusted for graphs
pp = None
adjust_cols_in_graph = ['n_frameshift_nested', 'n_strongBinders_nested']
values_to_plot = ['n_mutations_nested', 'n_missense_nested', 'n_frameshift_nested', 'n_binders_nested',
                  'n_weakBinders_nested', 'n_strongBinders_nested', 'n_summedBinders_nested', 'best_quality_nested',
                  'best_kDmt_nested']


# return integer range of numbers between min and max values in column
def integer_axis(a):
    return range(math.floor(min(graph_df[a])), math.ceil(max(graph_df[a]) + 1))


# save graph to pdf page
def save_to_pdf():
    print("entered save_to_pdf")
    pp.savefig(plt.gcf())
    plt.clf()


def draw_graph(nx_graph):
    fig, axes = plt.subplots(1,1,dpi=72)
    pos = graphviz_layout(nx_graph, prog="dot")
    nx.draw(nx_graph, pos, node_color=color_palette, ax=axes, with_labels=True)
    print("tree done")
    save_to_pdf()


def add_nodes(n, e, c):
    for index, row in graph_df.iterrows():
        n += [row['Clone_ID']]
        c += [float(row['CCF_X'])]
        if row['parent_clone_id'] != 0:
            e += [(row['parent_clone_id'], row['Clone_ID'])]
    return n, e, c


def draw_tree():
    global color_palette
    nodes = []
    edges = []
    color_gradient = []
    G = nx.Graph()
    nodes, edges, color_gradient = add_nodes(nodes, edges, color_gradient)
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["yellow", "orange", "crimson"])
    color_palette = [cmap(float(c/1.0)) for c in color_gradient]
    draw_graph(G)


# label clone points on the graph
def label_point(x_col, y_col, label, ax):
    a = pd.concat({'x': x_col, 'y': y_col, 'label': label}, axis=1)
    texts = []
    for i, point in a.iterrows():
        texts.append(ax.text(point['x'], point['y'], int(point['label'])))
    adjust_text(texts, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    save_to_pdf()


# draw line between parent and child point
def draw_line(c1, c2, x_axis, y_axis):
    x = []
    y = []
    new_df = graph_df[graph_df['Clone_ID'] == int(c1)]
    x.append(new_df[x_axis].iloc[0])
    y.append(new_df[y_axis].iloc[0])
    new_df = graph_df[graph_df['Clone_ID'] == int(c2)]
    x.append(new_df[x_axis].iloc[0])
    y.append(new_df[y_axis].iloc[0])
    plt.plot(x, y, color="b")


# graph using seaborn
def graph(x_axis, y_axis, var, sample):
    graph_df[x_axis] = [x for x in graph_df[x_axis] if x != 'nan']
    graph_df[y_axis] = [y for y in graph_df[y_axis] if y != 'nan']

    plot = sns.scatterplot(data=graph_df, x=x_axis, y=y_axis, hue=var, legend=False)
    plot.set_title(sample)

    for index, row in graph_df.iterrows():
        if row['parent_clone_id'] != 0:
            draw_line(row['parent_clone_id'], row['Clone_ID'], x_axis, y_axis)
    label_point(graph_df[x_axis], graph_df[y_axis], graph_df['Clone_ID'], plt.gca())


# plot vertical graph using CCF_X as the y-axis
def vertical_graph(sample):
    x_axis_list = values_to_plot
    y_axis = 'CCF_X'
    for x_axis in x_axis_list:
        print(x_axis)
        if x_axis in adjust_cols_in_graph:
            x_ticks = integer_axis(x_axis)
            plt.xticks(x_ticks)
        graph(x_axis, y_axis, 'Clone_ID', sample)


# plot horizontal graph using clonal level as the x-axis
def horizontal_graph(sample):
    x_axis = 'clonal_level'
    y_axis_list = values_to_plot
    for y_axis in y_axis_list:
        print(y_axis)
        x_ticks = integer_axis(x_axis)
        plt.xticks(x_ticks)
        if y_axis in adjust_cols_in_graph:
            y_ticks = integer_axis(y_axis)
            plt.yticks(y_ticks)
        graph(x_axis, y_axis, 'Clone_ID', sample)


# for each patient -- for each sample, create pdf object and plot graph
def start_graph():
    global clone_summary, graph_df, pp
    for p in patient_summary['Patient'].unique():
        print(p)
        for sample in clone_summary['Sample'].unique():
            if p in sample:
                pp = PdfPages(output_path + sample + '_figures.pdf')
                print(sample)
                graph_df = clone_summary[clone_summary['Sample'] == sample]
                draw_tree()
                vertical_graph(sample)
                horizontal_graph(sample)
                pp.close()


# for clone_summary file, add nested data columns for each clone
# important step for graphing clone data
def nested_data():
    global clone_summary
    clone_summary.best_quality.fillna(value=0, inplace=True)
    clone_summary.best_kDmt.fillna(value=100000, inplace=True)

    nested_cols = ['n_mutations_nested', 'n_missense_nested', 'n_frameshift_nested',
                                'n_binders_nested', 'n_weakBinders_nested', 'n_strongBinders_nested',
                                'n_summedBinders_nested', 'best_quality_nested',
                                'best_kDmt_nested', 'clonal_level']

    for index, row in clone_summary.iterrows():
        if row['parent_clone_id'] == 0:
            self_vals = [row['n_mutations'], row['n_missense'], row['n_frameshift'], row['n_neoantigen_binders'],
                         row['n_weakBinders'], row['n_strongBinders'], row['n_summedBinders'], row['best_quality'],
                         np.nan, 0]
            insert_value(clone_summary, index, nested_cols, self_vals)
        else:
            current_sample = row['Sample']
            parents = row['all_parent_clones'].split(",")[1:] + [str(row['Clone_ID'])]

            clonal_level = len(parents) - 1
            mutations_nested = 0
            missense_nested = 0
            frameshift_nested = 0
            binders_nested = 0
            weakBinders_nested = 0
            strongBinders_nested = 0
            quality_nested = 0
            kDmt_nested = 0

            for p in parents:
                temp_df = (clone_summary.loc[(clone_summary['Sample'] == current_sample) & (clone_summary['Clone_ID'] == int(p))]).iloc[0]

                mutations_nested += temp_df['n_mutations']
                missense_nested += temp_df['n_missense']
                frameshift_nested += temp_df['n_frameshift']
                binders_nested += temp_df['n_neoantigen_binders']
                weakBinders_nested += temp_df['n_weakBinders']
                strongBinders_nested += temp_df['n_strongBinders']

                if temp_df['best_quality']:
                    if ((temp_df['best_quality'] > quality_nested) | (quality_nested == 0)):
                        quality_nested = temp_df['best_quality']
                else:
                    quality_nested += 0
                if temp_df['best_kDmt']:
                    if ((temp_df['best_kDmt'] < kDmt_nested) | (kDmt_nested == 0)):
                        kDmt_nested = temp_df['best_kDmt']
                else:
                    kDmt_nested += 0

            nested_vals = [mutations_nested, missense_nested, frameshift_nested, binders_nested,
                           weakBinders_nested, strongBinders_nested, weakBinders_nested+strongBinders_nested,
                           quality_nested, kDmt_nested, clonal_level]
            insert_value(clone_summary, index, nested_cols, nested_vals)

    adjust_cols_in_clone = ['best_quality', 'best_kDmt', 'best_quality_nested', 'best_kDmt_nested']
    original_values = [0, 100000, 0, 100000]

    for i in range(0, len(adjust_cols_in_clone)):
        replace_in_col(clone_summary[adjust_cols_in_clone[i]], original_values[i], np.nan)
    clone_summary['n_summedBinders_nested'].replace(np.nan, 0, inplace=True)


# returns the path from root node to target node
def get_all_parents(json_tree, target_id):
    for element in json_tree:
        if 'clone_id' in element:
            if element['clone_id'] == target_id:
                return [element['clone_id']]
            else:
                if 'children' in element:
                    check_child = get_all_parents(element['children'], target_id)
                    if check_child:
                        return [element['clone_id']] + check_child


# search sample's phylogeny tree to collect mutation-specific and clone-specific data
# add data to mutation and clone summary files
def tree_data(tree, sample):
    global mutation_summary, clone_summary
    vcf_data = None
    neoantigen_sub_df = None
    for key in vcf_dict:
        if key in sample:
            vcf_data = vcf_dict[key]
            neoantigen_sub_df = neoantigen_dict[key]

    stack = tree['children'][:]
    while stack:
        item = stack.pop()

        parent_list = [0] + get_all_parents(tree['children'], item['clone_id'])

        missense_counter = 0
        frameshift_counter = 0
        mutations_counter = 0
        clone_neoantigens = 0
        clone_weakBinders = 0
        clone_strongBinders = 0

        for mutation_id in item["clone_mutations"]:
            mutations_counter += 1
            mutation_type = splice((vcf_data.loc[vcf_data.index[vcf_data['ID'] == mutation_id], 'ANN'].iloc[0])[0])

            if mutation_type.startswith("missense"):
                missense_counter += 1
            if mutation_type.startswith("frameshift"):
                frameshift_counter += 1

            sub_mutation_df = neoantigen_sub_df[neoantigen_sub_df['neoantigen'].str.startswith(mutation_id)]
            if len(sub_mutation_df) > 0:
                binders = len(sub_mutation_df)
                gene = sub_mutation_df['gene'].values[0]
                mutation_kdmt = sub_mutation_df.min(axis=0)['kDmt']
                mutation_strongBinders = len(sub_mutation_df[sub_mutation_df["kDmt"] < 50])
                mutation_weakBinders = len(sub_mutation_df[(sub_mutation_df["kDmt"] < 500)]) - mutation_strongBinders

                clone_neoantigens += binders
                clone_strongBinders += mutation_strongBinders
                clone_weakBinders += mutation_weakBinders

                mutation_summary = pd.concat([mutation_summary, pd.DataFrame({"Sample": [sample], "Mutation_ID": [mutation_id],
                                                                              "clone_id": item['clone_id'],
                                                                              "CCF_X": item['X'],
                                                                              "marginalCCF_x": item['x'],
                                                                              "mutation_type": mutation_type, "gene": gene,
                                                                              "n_neoantigen_binders": binders,
                                                                              "best_kDmt": mutation_kdmt,
                                                                              "n_weakBinders": mutation_weakBinders,
                                                                              "n_strongBinders": mutation_strongBinders})],
                                             ignore_index=True)

            else:
                mutation_summary = pd.concat([mutation_summary, pd.DataFrame({"Sample": [sample], "Mutation_ID": [mutation_id],
                                                                              "clone_id": item['clone_id'],
                                                                              "CCF_X": item['X'],
                                                                              "marginalCCF_x": item['x'],
                                                                              "mutation_type": mutation_type, "gene": None,
                                                                              "n_neoantigen_binders": 0,
                                                                              "best_kDmt": 0,
                                                                              "n_weakBinders": 0,
                                                                              "n_strongBinders": 0})],
                                             ignore_index=True)

        sub_clonal_df = filter_df_by_value(neoantigen_sub_df, 'clone_number', item['clone_id'], False)
        quality = sub_clonal_df.max(axis=0)['quality']
        fitness = sub_clonal_df.min(axis=0)['clone_fitness']
        clone_kdmt = sub_clonal_df.min(axis=0)['kDmt']
        clone_summary = pd.concat([clone_summary, pd.DataFrame({"Sample": [sample], "Clone_ID": [item['clone_id']],
                                                                "CCF_X": item['X'], "marginalCCF_x": item['x'],
                                                                "n_mutations": mutations_counter,
                                                                "n_missense": missense_counter,
                                                                "n_frameshift": frameshift_counter,
                                                                "parent_clone_id": parent_list[-2],
                                                                "all_parent_clones": ','.join(str(x) for x in parent_list[:-1]),
                                                                "n_neoantigen_binders": clone_neoantigens,
                                                                "n_weakBinders": clone_weakBinders,
                                                                "n_strongBinders": clone_strongBinders,
                                                                "n_summedBinders": clone_weakBinders+clone_strongBinders,
                                                                "clone_fitness": fitness,
                                                                "best_quality": quality,
                                                                "best_kDmt": clone_kdmt})],
                                  ignore_index=True)

        if 'children' in item:
            stack.extend(reversed(item['children']))
    return -1


# for each patient, open respective phylogeny tree json file
# save json tree in tree_dict with sample as key
# pass tree and sample name into tree_data function for data collection
def traverse_tree(name):
    tfile = "trees_" + name + ".json"
    f = open(tfile)
    data = json.load(f)

    original_tree = data['time_points'][0]['samples']
    for sub_tree in original_tree:
        print(sub_tree['id'])
        tree_dict[sub_tree['id']] = sub_tree['sample_trees'][0]['topology']
        tree_data(sub_tree['sample_trees'][0]['topology'], sub_tree['id'])
        print("\n")


# select patient tree and pass into traverse_tree function for sub_tree traversal
def tree_summary(path):
    os.chdir(path)
    print("entered")
    if patient:
        traverse_tree(patient)
    else:
        for p_name in patient_summary.Patient.unique():
            traverse_tree(p_name)


# add total mutations and total non-synonymous mutations columns to patient summary file
def total_cols():
    patient_summary.fillna(0, inplace=True)

    columns_list = list(patient_summary)
    col_list_non_syn = []

    non_syn = ["splice", "missense_variant", "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
               "disruptive_inframe_deletion", "disruptive_inframe_insertion"]
    for col in columns_list:
        if col.startswith(tuple(non_syn)):
            col_list_non_syn.append(col)
    patient_summary['non_syn_mutations'] = patient_summary[col_list_non_syn].sum(axis=1)

    columns_to_remove = ['Patient', 'Sample']
    for col in columns_to_remove:
        columns_list.remove(col)
    patient_summary['total_mutations'] = patient_summary[columns_list].sum(axis=1)


# save neo antigen dataframe in neoantigen_dict with sample as key
# count binders and collect data to add to patient summary file
def neoantigen_data(data, nf_filename):
    unique_9mers = []
    unique_binders = []
    patient_weakBinders = 0
    patient_strongBinders = 0

    for index, row in data.iterrows():
        if row['peptideMT'] not in unique_9mers:
            unique_9mers.append(row['peptideMT'])
        if row['neoantigen'] not in unique_binders:
            unique_binders.append(row['neoantigen'])
        if row['kDmt'] < 500:
            if row['kDmt'] > 50:
                patient_weakBinders += 1
            else:
                patient_strongBinders += 1

    for y in patient_summary['Sample']:
        if y in nf_filename:
            neoantigen_dict[y] = data
            index = patient_summary.index[patient_summary['Sample'] == y]
            nf_cols = ["n_unique_9mers", "n_neoantigen_binders", "n_weakBinders", "n_strongBinders"]
            nf_vals = [len(unique_9mers), len(unique_binders), patient_weakBinders, patient_strongBinders]
            insert_value(patient_summary, index, nf_cols, nf_vals)


# for each patient, open respective neo antigen file and read into dataframe
def neoantigen_summary(path):
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.startswith("nf_"):
                if patient:
                    if patient in f:
                        data = pd.read_csv(os.path.join(path, f), sep='\t', header=0)
                        neoantigen_data(data, f)
                else:
                    data = pd.read_csv(os.path.join(path, f), sep='\t', header=0)
                    neoantigen_data(data, f)
    total_cols()


# count unique mutations and feature types, add to patient summary file
def count(sample):
    info_list = vcf_dict[sample]["ANN"]

    feature_type = []
    for row in info_list:
        first_element = row[0]
        feature_type.append(splice(first_element))

    unique_list = []
    for mutation_type in feature_type:
        if mutation_type not in unique_list:
            unique_list.append(mutation_type)

    for unique_mutation in unique_list:
        index = patient_summary.index[patient_summary['Sample'] == sample]
        insert_value(patient_summary, index, [unique_mutation], [feature_type.count(unique_mutation)])


# read VCF files into pandas dataframe, store in vcf_dict with patient as key
def read_vcf(filepath, sample):
    reader = vcf.Reader(open(filepath))
    df = pd.DataFrame([vars(r) for r in reader])
    vcf_df = df.merge(pd.DataFrame(df.INFO.tolist()), left_index=True, right_index=True)
    vcf_dict[sample] = vcf_df
    count(sample)


# within patient VCF folder, search for VCF files
def search_dir(path, p):
    global patient_summary
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.endswith(vcf_suffix):
                sample_name = f[:f.find(vcf_suffix)]
                patient_summary = pd.concat([patient_summary, pd.DataFrame({"Patient": [p], "Sample": [sample_name]})],
                                            ignore_index=True)
                read_vcf(os.path.join(path, f), sample_name)


# given root directory, find VCF directory and loop through patient folders
def findpath(path):
    os.chdir(path)
    if os.path.isdir("VCF"):
        os.chdir(path + "/VCF")
        if patient:
            new_path = os.path.join(os.getcwd(), patient)
            search_dir(new_path, patient)
        else:
            for sub_dir in os.listdir():
                new_path = os.path.join(os.getcwd(), sub_dir)
                search_dir(new_path, sub_dir)
    else:
        return "VCF not found in directory"


print("start")
# call vcf function, neo antigen function, phylogeny tree function
findpath(vcf_directory)
neoantigen_summary(neoant_directory)
tree_summary(tree_directory)

# If User specified patient, double check that output files ONLY contain patient-specific data
if patient:
    patient_summary = filter_df_by_value(patient_summary, 'Patient', patient, False)
    mutation_summary = filter_df_by_value(mutation_summary, 'Sample', patient, True)
    clone_summary = filter_df_by_value(clone_summary, 'Sample', patient, True)

# replace 0 values with N/A in best_kDmt column
replace_in_col(mutation_summary['best_kDmt'], 0, np.nan)

# call function to add nested columns to clone_summary file
nested_data()

# call graph function to draw seaborn plots and phylogeny tree
start_graph()

# export summary files as .csv
patient_summary.to_csv(output_path + "082224_patient_summary.csv")
clone_summary.to_csv(output_path + "082224_clone_summary.csv")
mutation_summary.to_csv(output_path + "082224_mutation_summary.csv")

print("end")
