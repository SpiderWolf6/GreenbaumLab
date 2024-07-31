import io, os
import pandas as pd
import vcf, vcf.parser
import sys
import requests, json, argparse

parser = argparse.ArgumentParser(description='Foo')
parser.add_argument('--vcf_dir', metavar='vcf_dir', type=str, required=True)
parser.add_argument('--vcf_suff', metavar='vcf_suff', type=str, required=True)
parser.add_argument('--tree_dir', metavar='tree_dir', type=str, required=True)
parser.add_argument('--na_dir', metavar='na_dir', type=str, required=True)
parser.add_argument('--patient', metavar='patient', type=str, required=False)
args = parser.parse_args()

vcf_directory = args.vcf_dir
patient = args.patient
vcf_suffix = args.vcf_suff
neoant_directory = args.na_dir
tree_directory = args.tree_dir

patient_summary = pd.DataFrame(columns=['Patient', 'Sample'], index=[])
mutation_summary = pd.DataFrame()
clone_summary = pd.DataFrame()

vcf_dict = {}
neoantigen_dict = {}
tree_dict = {}

def traverse_tree(tree, sample):
    global mutation_summary
    global clone_summary
    vcf_data = None
    neoantigen_data = None
    for key in vcf_dict:
        if key in sample:
            vcf_data = vcf_dict[key]
            neoantigen_data = neoantigen_dict[key]

    stack = tree['children'][:]
    parents = []
    all_parents = []
    while stack:
        parents.clear()
        all_parents.clear()
        item = stack.pop()

        missense_counter = 0
        frameshift_counter = 0
        mutations_counter = 0

        for m in item["clone_mutations"]:
            mutations_counter += 1
            mut_type = (vcf_data.loc[vcf_data.index[vcf_data['ID'] == m], 'ANN'].iloc[0])[0]
            start = mut_type.find("|") + 1
            end = mut_type[start:].find("|") + start
            mut_type = mut_type[start:end]

            if mut_type.startswith("missense"):
                missense_counter += 1
            if mut_type.startswith("frameshift"):
                frameshift_counter += 1

            sub_df = neoantigen_data[neoantigen_data['neoantigen'].str.startswith(m)]
            if len(sub_df) > 0:
                binders = len(sub_df)
                gene = sub_df['gene'].values[0]
                kDmt_min = sub_df.min(axis=0)['kDmt']

                mutation_summary = pd.concat([mutation_summary, pd.DataFrame({"Sample": [sample], "Mutation_ID": [m],
                                             "clone_id": item['clone_id'], "CCF_X": item['X'], "marginalCCF_x": item['x'],
                                             "mutation_type": mut_type, "gene": gene, "n_neoantigen_binders": binders,
                                             "best_kDmt": kDmt_min})],
                                             ignore_index=True)

            else:
                mutation_summary = pd.concat([mutation_summary, pd.DataFrame({"Sample": [sample], "Mutation_ID": [m],
                                             "clone_id": item['clone_id'], "CCF_X": item['X'], "marginalCCF_x": item['x'],
                                             "mutation_type": mut_type, "gene": None, "n_neoantigen_binders": 0,
                                             "best_kDmt": 0})],
                                             ignore_index=True)

        clone_summary = pd.concat([clone_summary, pd.DataFrame({"Sample": [sample], "Clone_ID": [item['clone_id']],
                                                                "CCF_X": item['X'], "marginalCCF_x": item['x'],
                                                                "n_mutations": mutations_counter,
                                                                "n_missense": missense_counter,
                                                                "n_frameshift": frameshift_counter,
                                                                "parent_clone_id": None, "all_parent_clones": None})],
                                                                ignore_index=True)

        if 'children' in item:
            stack.extend(reversed(item['children']))
    return -1

def tree_summary(path):
    os.chdir(path)
    print("entered")
    for p_name in patient_summary.Patient.unique():
        tfile = "trees_" + p_name + ".json"
        f = open(tfile)
        data = json.load(f)
        original_tree = data['time_points'][0]['samples']
        for x in original_tree:
            print(x['id'])
            tree_dict[x['id']] = x['sample_trees'][0]['topology']
            traverse_tree(x['sample_trees'][0]['topology'], x['id'])
            print("\n")

def total_cols():
    patient_summary.fillna(0, inplace=True)

    col_list = list(patient_summary)
    col_list_non_syn = []

    non_syn = ["splice", "missense_variant", "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
               "disruptive_inframe_deletion", "disruptive_inframe_insertion"]
    for x in col_list:
        if x.startswith(tuple(non_syn)):
            col_list_non_syn.append(x)
    patient_summary['non_syn_mutations'] = patient_summary[col_list_non_syn].sum(axis=1)

    columns_to_remove = ['Patient', 'Sample']
    for x in columns_to_remove:
        col_list.remove(x)
    patient_summary['total_mutations'] = patient_summary[col_list].sum(axis=1)

def neoantigen_summary(path):
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.startswith("nf_"):
                unique_9mers = []
                unique_binders = []
                weakbind_count = 0
                strongbind_count = 0

                data = pd.read_csv(os.path.join(path, f), sep='\t', header=0)
                for index, row in data.iterrows():
                    if row['peptideMT'] not in unique_9mers:
                        unique_9mers.append(row['peptideMT'])
                    if row['neoantigen'] not in unique_binders:
                        unique_binders.append(row['neoantigen'])
                    if row['kDmt'] < 500:
                        if row['kDmt'] > 50:
                            weakbind_count += 1
                        else:
                            strongbind_count += 1

                for y in patient_summary['Sample']:
                    if y in f:
                        neoantigen_dict[y] = data
                        i = patient_summary.index[patient_summary['Sample'] == y]
                        patient_summary.loc[i,["n_unique_9mers", "n_binders_9mers", "n_weakBinders_9mer", "n_strongBinders_9mer"]] = [len(unique_9mers), len(unique_binders), weakbind_count, strongbind_count]
    total_cols()

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

def read_vcf(filepath, sample):
    reader = vcf.Reader(open(filepath))
    df = pd.DataFrame([vars(r) for r in reader])
    vcf_df = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    vcf_dict[sample] = vcf_df
    count(sample)

# within
def search_dir(path, p):
    global patient_summary
    for dirpath, dirname, filename in os.walk(path):
        for f in filename:
            if f.endswith(vcf_suffix):
                sample_name = f[:f.find(vcf_suffix)]
                patient_summary = pd.concat([patient_summary, pd.DataFrame({"Patient": [p], "Sample": [sample_name]})], ignore_index=True)
                read_vcf(os.path.join(path, f), sample_name)

# given root directory, find VCF directory or specific patient directory
def findpath(path, p):
    os.chdir(path)
    if os.path.isdir("VCF"):
        os.chdir(path + "/VCF")
        if p:
            print("patient is given")
            if os.path.isdir(p):
                os.chdir(os.getcwd() + "/" + p)
                search_dir(os.getcwd(), p)
        else:
            for sub_dir in os.listdir():
                new_path = os.path.join(os.getcwd(), sub_dir)
                search_dir(new_path, sub_dir)
    else:
        return "VCF not found in directory"

findpath(vcf_directory, patient)
neoantigen_summary(neoant_directory)
print(patient_summary)


tree_summary(tree_path)
mutation_summary.to_csv()
patient_summary.to_csv()
