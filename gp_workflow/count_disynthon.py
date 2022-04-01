import numpy as np
import pandas as pd
import pickle
import argparse

def count_disynthon(comp, bb_dict):
    # Load in dictionary to convert building block SMILES names into hashes
    # to make it easier to keep track of

    compounds = pd.DataFrame()
    compounds['bb1'] = comp['bb1'].apply(lambda x: bb_dict[x])
    compounds['bb2'] = comp['bb2'].apply(lambda x: bb_dict[x])
    compounds['bb3'] = comp['bb3'].apply(lambda x: bb_dict[x])
    compounds['active'] = comp['active']

    ## List out the disynthon pairs as tuple so they can be indexed in a dictionary
    compounds['AB'] = [(x,y) for x,y in zip(compounds['bb1'], compounds['bb2'])]
    compounds['BC'] = [(y,z) for y,z in zip(compounds['bb2'], compounds['bb3'])]
    compounds['AC'] = [(x,z) for x,z in zip(compounds['bb1'], compounds['bb3'])]

    active = compounds.loc[compounds['active'] == True]
    inactive = compounds.loc[compounds['active'] == False]

    ## Count how many times each disynthon occurs in actives and inactives
    disynthon_dict = {}

    ## Keeping track of disynthons in active compounds
    for disynthon in active['AB']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [1, 0]
        else:
            disynthon_dict[disynthon][0] += 1

    for disynthon in active['BC']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [1, 0]
        else:
            disynthon_dict[disynthon][0] += 1

    for disynthon in active['AC']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [1, 0]
        else:
            disynthon_dict[disynthon][0] += 1

    ## Keeping track of disynthons in inactive compounds
    for disynthon in inactive['AB']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [0, 1]
        else:
            disynthon_dict[disynthon][1] += 1

    for disynthon in inactive['BC']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [0, 1]
        else:
            disynthon_dict[disynthon][1] += 1

    for disynthon in inactive['AC']:
        if disynthon not in disynthon_dict:
            disynthon_dict[disynthon] = [0, 1]
        else:
            disynthon_dict[disynthon][1] += 1

    return disynthon_dict


if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(
    description="Get counts of how many times disynthons appear in actives and inactives",
    allow_abbrev=False)

    my_parser.add_argument('--compounds',
            action='store',
            type=str,
            help='path to pandas dataframe of compounds',
            required=True)
    my_parser.add_argument('--dict',
            action='store',
            type=str,
            help='path to dictionary of building block hashes',
            required=True)

    args = my_parser.parse_args()

    disynth_dict = count_disynthon(args.compounds, args.dict)
    disynth_dict.to_csv('{}_disynth.csv'.format(args.compounds[:-4]), index=False)
