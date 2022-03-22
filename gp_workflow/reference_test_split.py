import numpy as np
import pandas as pd
from sys import argv
import argparse
from count_disynthon import *

def ref_test_split(seed, path):
    '''
    Creates a random split of the angx sEH dataset based on the random seed provided.
    '''
    path_to_actives = '../data/del_actives_no_boron.csv'
    path_to_inactives = '../data/del_inactives_no_boron.csv'
    path_to_bbs = '../data/deprotected_bbs_updated.csv'
    path_to_dict = '../data/deprotected_bbs_dictionary_updated.pkl'

    ## Create a concatenated table of all compounds from the DEL
    actives = pd.read_csv(path_to_actives).drop(columns=['read_count'])
    actives['active'] = [True for x in range(len(actives))]
    inactives = pd.read_csv(path_to_inactives)
    inactives['active'] = [False for x in range(len(inactives))]
    all_compounds = pd.concat([actives, inactives], ignore_index=True)

    ## Choose how many building blocks to look across for a reference and test split
    p = 0.6
    N = 1500
    full_bb_directory = pd.read_csv(path_to_bbs)
    full_bb_dictionary = pickle.load(open(path_to_dict, 'rb'))
    selected_bbs = full_bb_directory[:N]
    ref_frac = int(p*len(selected_bbs))

    ## Set the random seed based on input flag
    np.random.seed(seed)
    split = np.random.permutation(selected_bbs['bb'])
    reference_bbs = split[:ref_frac]
    test_bbs = split[ref_frac:]

    ##Query all compounds to see which of its building blocks belong in either the
    ##defined reference or test set -- throw out compounds with some building blocks
    ##in one and the rest in the other
    all_compounds['bb1_ref'] = all_compounds['bb1'].apply(lambda x: x in reference_bbs)
    all_compounds['bb2_ref'] = all_compounds['bb2'].apply(lambda x: x in reference_bbs)
    all_compounds['bb3_ref'] = all_compounds['bb3'].apply(lambda x: x in reference_bbs)

    all_compounds['bb1_test'] = all_compounds['bb1'].apply(lambda x: x in test_bbs)
    all_compounds['bb2_test'] = all_compounds['bb2'].apply(lambda x: x in test_bbs)
    all_compounds['bb3_test'] = all_compounds['bb3'].apply(lambda x: x in test_bbs)

    ref_index = list(set(np.where(all_compounds[['bb1_ref', 'bb2_ref', 'bb3_ref']].sum(axis=1) == 3)[0]) & set(np.where(all_compounds[['bb1_test', 'bb2_test', 'bb3_test']].sum(axis=1) == 0)[0]))
    test_index = list(set(np.where(all_compounds[['bb1_test', 'bb2_test', 'bb3_test']].sum(axis=1) == 3)[0]) & set(np.where(all_compounds[['bb1_ref', 'bb2_ref', 'bb3_ref']].sum(axis=1) == 0)[0]))

    ref_compounds = all_compounds.iloc[ref_index]
    test_compounds = all_compounds.iloc[test_index]

    ## Perform disynthon aggregation and return disynthons and their counts
    ref_disynth_dict = count_disynthon(ref_compounds, full_bb_dictionary)
    pickle.dump(ref_disynth_dict, open('{}/reference_disynth_dict_{}.pkl'.format(path, seed), 'wb'))

    test_disynth_dict = count_disynthon(test_compounds, full_bb_dictionary)
    pickle.dump(test_disynth_dict, open('{}/test_disynth_dict_{}.pkl'.format(path, seed), 'wb'))

    ## Save the reference and test compound SMILES
    ref_compounds['SMILES'].to_csv('{}/reference_compounds_{}.csv'.format(path, seed), index=False)
    test_compounds['SMILES'].to_csv('{}/test_compounds_{}.csv'.format(path, seed), index=False)

    ## Extract all building blocks from reference and test set splits
    ref_bbs = pd.DataFrame(columns=['SMILES'])
    ref_bbs['SMILES'] = list(set(ref_compounds['bb1']) | set(ref_compounds['bb2']) | set(ref_compounds['bb3']))
    ref_bbs.to_csv('{}/reference_bbs_{}.csv'.format(path, seed), index=False)

    test_bbs = pd.DataFrame(columns=['SMILES'])
    test_bbs['SMILES'] = list(set(test_compounds['bb1']) | set(test_compounds['bb2']) | set(test_compounds['bb3']))
    test_bbs.to_csv('{}/test_bbs_{}.csv'.format(path, seed), index=False)


if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(
    description="Create reference and test set split from list of compounds",
    allow_abbrev=False)

    my_parser.add_argument('--seed',
            action='store',
            type=int,
            help='value of random seed',
            required=True)

    my_parser.add_argument('--path',
            action='store',
            type=str,
            help='path to save files',
            required=True)

    args = my_parser.parse_args()

    ref_test_split(args.seed, args.path)
