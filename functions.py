import os
import math
import json
import vcf
import pydot

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
from adjustText import adjust_text

from functools import wraps

# for seaborn plots, change size of canvas and create color palette object
sns.set_theme(rc={'figure.figsize': (11.7, 8.27)})
color_palette = None

os.environ["PATH"] += r"" # local path to Graphviz/bin

# create vcf, neoantigen, and tree dictionaries to organize individual datasets
vcf_dict = {}
neoantigen_dict = {}
tree_dict = {}

# create pdf object, and the axes that will need to plotted and adjusted for graphs
pp = None
adjust_cols_in_graph = ['n_frameshift_nested', 'n_strongBinders_nested']
values_to_plot = ['n_mutations_nested', 'n_missense_nested', 'n_frameshift_nested', 'n_binders_nested',
                  'n_weakBinders_nested', 'n_strongBinders_nested', 'n_summedBinders_nested', 'best_quality_nested',
                  'best_kDmt_nested']

def debug(func):
    wraps(func)
    def wrapper(*args, **kwargs):
        print(f"Calling {func.__name__} with args: {args} and kwargs: {kwargs}")
        result = func(*args, **kwargs)
        print(f"{func.__name__} returned: {result}")
        return result
    return wrapper

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

# return integer range of numbers between min and max values in column
def integer_axis(df, a):
    return range(math.floor(min(df[a])), math.ceil(max(df[a]) + 1))

# save graph to pdf page
def save_to_pdf():
    print("entered save_to_pdf")
    pp.savefig(plt.gcf())
    plt.clf()

def draw_graph(nx_graph):
    print("call draw_graph, draw tree diagram")
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
    print("call draw_tree, create nodes and edges")
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
    print("call graph, draw graphs for " + sample)
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
            x_ticks = integer_axis(graph_df, x_axis)
            plt.xticks(x_ticks)
            graph(x_axis, y_axis, 'Clone_ID', sample)

# plot horizontal graph using clonal level as the x-axis
def horizontal_graph(sample):
    x_axis = 'clonal_level'
    y_axis_list = values_to_plot
    for y_axis in y_axis_list:
        print(y_axis)
        x_ticks = integer_axis(graph_df, x_axis)
        plt.xticks(x_ticks)
        if y_axis in adjust_cols_in_graph:
            y_ticks = integer_axis(graph_df, y_axis)
            plt.yticks(y_ticks)
        graph(x_axis, y_axis, 'Clone_ID', sample)

# for each patient -- for each sample, create pdf object and plot graph
def start_graph(patient_df, clone_df, output_path):
    print("call start_graph")
    global graph_df, pp
    for p in patient_df['Patient'].unique():
        print(p)
        for sample in clone_df['Sample'].unique():
            if p in sample:
                pp = PdfPages(output_path + '/' + sample + '_figures.pdf')
                print(sample)
                graph_df = clone_df[clone_df['Sample'] == sample]
                draw_tree()
                vertical_graph(sample)
                horizontal_graph(sample)
                pp.close()
    return "drawing trees done"

# for clone_summary file, add nested data columns for each clone
# important step for graphing clone data
def nested_data(df):
    print("call nested_data, collecting child clone data")
    df.best_quality.fillna(value=0, inplace=True)
    df.best_kDmt.fillna(value=100000, inplace=True)

    nested_cols = ['n_mutations_nested', 'n_missense_nested', 'n_frameshift_nested',
                                'n_binders_nested', 'n_weakBinders_nested', 'n_strongBinders_nested',
                                'n_summedBinders_nested', 'best_quality_nested',
                                'best_kDmt_nested', 'clonal_level']

    for index, row in df.iterrows():
        if row['parent_clone_id'] == 0:
            self_vals = [row['n_mutations'], row['n_missense'], row['n_frameshift'], row['n_neoantigen_binders'],
                         row['n_weakBinders'], row['n_strongBinders'], row['n_summedBinders'], row['best_quality'],
                         np.nan, 0]
            insert_value(df, index, nested_cols, self_vals)
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
                temp_df = (df.loc[(df['Sample'] == current_sample) & (df['Clone_ID'] == int(p))]).iloc[0]

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
            insert_value(df, index, nested_cols, nested_vals)

    adjust_cols_in_clone = ['best_quality', 'best_kDmt', 'best_quality_nested', 'best_kDmt_nested']
    original_values = [0, 100000, 0, 100000]

    for i in range(0, len(adjust_cols_in_clone)):
        replace_in_col(df[adjust_cols_in_clone[i]], original_values[i], np.nan)
    df['n_summedBinders_nested'].replace(np.nan, 0, inplace=True)
    return df

# returns the path from root node to target node
def get_all_parents(json_tree, target_id):
    print("call get_all_parents, retrieve path to node")
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
def tree_data(tree, sample, mutation_df, clone_df):
    print("call tree_data, traversing sub_trees for " + sample)
    vcf_data = None
    neoantigen_sub_df = None
    for key in vcf_dict:
        print(key)
        print(sample)
        print("neoantigen dict time")
        if key in sample:

            vcf_data = vcf_dict[key]
            print(neoantigen_dict)
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

                mutation_df = pd.concat([mutation_df, pd.DataFrame({"Sample": [sample], "Mutation_ID": [mutation_id],
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
                mutation_df = pd.concat([mutation_df, pd.DataFrame({"Sample": [sample], "Mutation_ID": [mutation_id],
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
        clone_df = pd.concat([clone_df, pd.DataFrame({"Sample": [sample], "Clone_ID": [item['clone_id']],
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
    print("burger")
    print(clone_df)
    return mutation_df, clone_df

# for each patient, open respective phylogeny tree json file
# save json tree in tree_dict with sample as key
# pass tree and sample name into tree_data function for data collection
def traverse_tree(patient):
    print("call traverse_tree, reading " + patient + " tree")
    tfile = "trees_" + patient + ".json"
    f = open(tfile)
    data = json.load(f)

    original_tree = data['time_points'][0]['samples']
    for sub_tree in original_tree:
        print(sub_tree['id'])
        tree_dict[sub_tree['id']] = sub_tree['sample_trees'][0]['topology']
        child_tree = sub_tree['sample_trees'][0]['topology']
        tree_id = sub_tree['id']
        return child_tree, tree_id

# select patient tree and pass into traverse_tree function for sub_tree traversal
def tree_summary(patient_df, path, given_patient):
    print("call tree_summary, entered trees folder")
    os.chdir(path)
    if given_patient:
        return given_patient
    else:
        for patient in patient_df.Patient.unique():
            return patient

# add total mutations and total non-synonymous mutations columns to patient summary file
def total_cols(df):
    print("call total_cols, creating total columns for patient summary file")
    df.fillna(0, inplace=True)

    columns_list = list(df)
    col_list_non_syn = []

    non_syn = ["splice", "missense_variant", "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
               "disruptive_inframe_deletion", "disruptive_inframe_insertion"]
    for col in columns_list:
        if col.startswith(tuple(non_syn)):
            col_list_non_syn.append(col)
    df['non_syn_mutations'] = df[col_list_non_syn].sum(axis=1)

    columns_to_remove = ['Patient', 'Sample']
    for col in columns_to_remove:
        columns_list.remove(col)
    df['total_mutations'] = df[columns_list].sum(axis=1)
    return df

# save neo antigen dataframe in neoantigen_dict with sample as key
# count binders and collect data to add to patient summary file
def neoantigen_data(patient_df, df, filename):
    print("call neoantigen_data, collecting neoantigen data for " + filename)
    unique_9mers = []
    unique_binders = []
    patient_weakBinders = 0
    patient_strongBinders = 0

    for index, row in df.iterrows():
        if row['peptideMT'] not in unique_9mers:
            unique_9mers.append(row['peptideMT'])
        if row['neoantigen'] not in unique_binders:
            unique_binders.append(row['neoantigen'])
        if row['kDmt'] < 500:
            if row['kDmt'] > 50:
                patient_weakBinders += 1
            else:
                patient_strongBinders += 1

    for y in patient_df['Sample']:
        if y in filename:

            print("the y sample is" + y)
            neoantigen_dict[y] = df
            index = patient_df.index[patient_df['Sample'] == y]
            nf_cols = ["n_unique_9mers", "n_neoantigen_binders", "n_weakBinders", "n_strongBinders"]
            nf_vals = [len(unique_9mers), len(unique_binders), patient_weakBinders, patient_strongBinders]
            insert_value(patient_df, index, nf_cols, nf_vals)
    return patient_df

# for each patient, open respective neo antigen file and read into dataframe
def neoantigen_summary(patient_df, path, given_patient):
    new_df = pd.DataFrame()
    print("call neoantigen_summary, entered neoantigen folder")
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.startswith("nf_"):
                if given_patient:
                    if given_patient in f:
                        df = pd.read_csv(os.path.join(path, f), sep='\t', header=0)
                        new_df = neoantigen_data(patient_df, df, f)
                else:
                    df = pd.read_csv(os.path.join(path, f), sep='\t', header=0)
                    new_df = neoantigen_data(patient_df, df, f)
    return new_df
# count unique mutations and feature types, add to patient summary file
def count(df, sample):
    print("call count, counting mutations for " + sample)
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
        index = df.index[df['Sample'] == sample]
        insert_value(df, index, [unique_mutation], [feature_type.count(unique_mutation)])
    return df

# read VCF files into pandas dataframe, store in vcf_dict with patient as key
def read_vcf(filepath, sample):
    global vcf_dict
    print("call read_vcf, entered " + sample + " folder")
    reader = vcf.Reader(open(filepath))
    df = pd.DataFrame([vars(r) for r in reader])
    vcf_df = df.merge(pd.DataFrame(df.INFO.tolist()), left_index=True, right_index=True)
    vcf_dict[sample] = vcf_df
    return sample

all_filepaths = []
all_samples = []

# within patient VCF folder, search for VCF files
def search_dir(df, path, patient, endString):
    print("call search_dir, entered 'VCF' folder, loop through patient folders")
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.endswith(endString):
                sample_name = f[:f.find(endString)]
                df = pd.concat([df, pd.DataFrame({"Patient": [patient], "Sample": [sample_name]})],
                                            ignore_index=True)
                all_filepaths.append(os.path.join(path,f))
                all_samples.append(sample_name)
                return df

# given root directory, find VCF directory and loop through patient folders
def findpath(df, path, given_patient,suffix):
    print("call findpath, searching for 'VCF' folder")
    os.chdir(path)
    if os.path.isdir("VCF"):
        print("'VCF' folder found")
        os.chdir(path + "/VCF")
        if given_patient:
            new_path = os.path.join(os.getcwd(), given_patient)
            df = search_dir(df, new_path, given_patient, suffix)
        else:
            for sub_dir in os.listdir():
                new_path = os.path.join(os.getcwd(), sub_dir)
                df = search_dir(df, new_path, sub_dir, suffix)
        return df, all_filepaths, all_samples
    else:
        return "VCF not found in directory"
