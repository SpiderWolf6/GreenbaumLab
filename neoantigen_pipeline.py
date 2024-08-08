import io
import os
import sys
import time

import pandas as pd
import vcf
import vcf.parser
import argparse
import requests
import json
import numpy as np

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
    neoantigen_data = None
    for key in vcf_dict:
        if key in sample:
            vcf_data = vcf_dict[key]
            neoantigen_data = neoantigen_dict[key]

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
            mutation_type = (vcf_data.loc[vcf_data.index[vcf_data['ID'] == mutation_id], 'ANN'].iloc[0])[0]
            start = mutation_type.find("|") + 1
            end = mutation_type[start:].find("|") + start
            mutation_type = mutation_type[start:end]

            if mutation_type.startswith("missense"):
                missense_counter += 1
            if mutation_type.startswith("frameshift"):
                frameshift_counter += 1

            sub_mutation_df = neoantigen_data[neoantigen_data['neoantigen'].str.startswith(mutation_id)]
            if len(sub_mutation_df) > 0:
                binders = len(sub_mutation_df)
                gene = sub_mutation_df['gene'].values[0]
                mutation_kDmt = sub_mutation_df.min(axis=0)['kDmt']
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
                                                                              "best_kDmt": mutation_kDmt,
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

        sub_clonal_df = neoantigen_data[neoantigen_data['clone_number'] == item['clone_id']]
        quality = sub_clonal_df.max(axis=0)['quality']
        fitness = sub_clonal_df.min(axis=0)['clone_fitness']
        clone_kDmt = sub_clonal_df.min(axis=0)['kDmt']
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
                                                                "clone_fitness": fitness,
                                                                "best_quality": quality,
                                                                "best_kDmt": clone_kDmt})],
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
        columns_list_list.remove(col)
    patient_summary['total_mutations'] = patient_summary[col_list].sum(axis=1)


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
            i = patient_summary.index[patient_summary['Sample'] == y]
            patient_summary.loc[i, ["n_unique_9mers", "n_neoantigen_binders", "n_weakBinders",
                                    "n_strongBinders"]] = [len(unique_9mers), len(unique_binders),
                                                                patient_weakBinders, patient_strongBinders]


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
        row = row[0]
        start = row.find("|") + 1
        end = row[start:].find("|") + start
        feature_type.append(row[start:end])

    unique_list = []

    for obj in feature_type:
        if obj not in unique_list:
            unique_list.append(obj)

    for obj in unique_list:
        patient_summary.loc[patient_summary.index[patient_summary['Sample'] == sample], obj] = feature_type.count(obj)


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


# call vcf function, neo antigen function, phylogeny tree function
findpath(vcf_directory)
neoantigen_summary(neoant_directory)
tree_summary(tree_directory)

# If User specified patient, check that summary files only contain patient-specific data
if patient:
    patient_summary = patient_summary[patient_summary['Patient'] == patient]
    mutation_summary = mutation_summary[mutation_summary['Sample'].str.contains(patient)]
    clone_summary = clone_summary[clone_summary['Sample'].str.contains(patient)]

mutation_summary['best_kDmt'].replace(0, np.nan, inplace=True)

mutation_summary.to_csv()
clone_summary.to_csv()
patient_summary.to_csv()
# export summary files as .csv
