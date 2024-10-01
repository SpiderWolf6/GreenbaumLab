import pandas as pd
import requests, argparse
from functions import *

os.environ["PATH"] += r"" # local path to Graphviz/bin

# initialize user input parser and arguments
parser = argparse.ArgumentParser(description='Foo')
parser.add_argument('--vcf_dir', metavar='vcf_dir', type=str, required=True)
parser.add_argument('--vcf_suff', metavar='vcf_suff', type=str, required=True)
parser.add_argument('--tree_dir', metavar='tree_dir', type=str, required=True)
parser.add_argument('--neoag_dir', metavar='neoag_dir', type=str, required=True)
parser.add_argument('--output_dir', metavar='output_dir', type=str, required=True)
parser.add_argument('--patient', metavar='patient', type=str, required=False)
args = parser.parse_args()

# convert arguments into global variables
vcf_directory = args.vcf_dir
patient = args.patient
vcf_suffix = args.vcf_suff
neoantigen_directory = args.neoag_dir
tree_directory = args.tree_dir
output_directory = args.output_dir

# create patient, mutation, and clone summary dataframes
patient_summary = pd.DataFrame(columns=['Patient', 'Sample'], index=[])
mutation_summary = pd.DataFrame()
clone_summary = pd.DataFrame()

# graph_df object for storing temporary dataframes used in plotting graphs and drawing trees
graph_df = None

# output folder
output_path = output_directory

files = []
samples = []

def main():
    global patient_summary, mutation_summary, clone_summary, files, samples
    print("start")
    # call vcf function, neo antigen function, phylogeny tree function
    patient_summary, files, samples = findpath(patient_summary, vcf_directory, patient, vcf_suffix)
    for index, e in enumerate(files):
        sample_name = read_vcf(files[index], samples[index])
        patient_summary = count(patient_summary, sample_name)
        print("vcf data analysis done")

    patient_summary = neoantigen_summary(patient_summary, neoantigen_directory, patient)
    print("neoantigen data analysis done")
    print(patient_summary)
    patient_summary = total_cols(patient_summary)
    print("total_cols added, patient summary table created")
    for p in patient_summary.Patient.unique():
        tree_patient = tree_summary(patient_summary, tree_directory, p)
        sub_tree, tree_id = traverse_tree(tree_patient)
        mutation_summary, clone_summary = tree_data(sub_tree, tree_id, mutation_summary, clone_summary)
    print("clone and mutation data added, all three summary files created")
    # print(patient_summary, clone_summary, mutation_summary)

    # If User specified patient, double check that output files ONLY contain patient-specific data
    if patient:
        patient_summary = filter_df_by_value(patient_summary, 'Patient', patient, False)
        mutation_summary = filter_df_by_value(mutation_summary, 'Sample', patient, True)
        clone_summary = filter_df_by_value(clone_summary, 'Sample', patient, True)

    # replace 0 values with N/A in best_kDmt column
    replace_in_col(mutation_summary['best_kDmt'], 0, np.nan)

    clone_summary = nested_data(clone_summary)
    print("nested data added to clone_file")

    graph_return_statement = start_graph(patient_summary, clone_summary, output_path)
    print(graph_return_statement)

    # export summary files as .csv
    patient_summary.to_csv(output_path + "/100124_patient_summary.csv")
    clone_summary.to_csv(output_path + "/100124_clone_summary.csv")
    mutation_summary.to_csv(output_path + "/100124_mutation_summary.csv")

    print("end")

main()
