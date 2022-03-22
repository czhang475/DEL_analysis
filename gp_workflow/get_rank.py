import numpy as np
import pandas as pd
import argparse

def get_rank(k, sim_matrix, test_compounds, output):
    sim_matrix = np.load(sim_matrix)
    test_compounds = pd.read_csv(test_compounds)
    # Sum the top k scores for each column
    top_k_scores = np.sort(sim_matrix, axis=0)[::-1][:k, :]
    sum_top = np.sum(top_k_scores, axis=0)
    # Rank column sums in descending order
    ranking = np.argsort(sum_top)[::-1]
    # Plug rankings back into the test dataframe to get ranks for compounds
    test_compounds.iloc[ranking].to_csv(output, index=False)

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(
    description="Extracts the top k scores of each column in a similarity matrix to rank compounds by decreasing likelihood of activity",
    allow_abbrev=False)

    my_parser.add_argument('--k',
            action='store',
            type=int,
            help='how many of the top compounds to use for scoring',
            required=True)

    my_parser.add_argument('--matrix',
            action='store',
            type=str,
            help='path to 3D similarity matrix',
            required=True)

    my_parser.add_argument('--test',
            action='store',
            type=str,
            help='path to .csv of test compounds',
            required=True)

    my_parser.add_argument('--output',
            action='store',
            type=str,
            help='output path of test compounds ranked',
            required=True)

    args = my_parser.parse_args()

    compound_rank = get_rank(args.k, args.matrix, args.test, args.output)
