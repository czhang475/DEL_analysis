from openeye import oechem, oefastrocs
import pickle
import os
import sys
import argparse
import numpy as np
import time
from clean_3D_sim_matrix import *

def timer(func):
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        statement = f'{total_time:.3f} seconds'
        with open('outfile/sim_3D_times.txt', 'a+') as f:
            f.write('\n' + statement)
        return value
    return wrapper_timer

def calc_3D_sim(ref, test, ref_groups, test_groups, output):
    if not oefastrocs.OEFastROCSIsGPUReady():
        oechem.OEThrow.Info("No supported GPU available!")

    opts = oefastrocs.OEShapeDatabaseOptions()
    print(opts.GetLimit())
    print(opts.GetScoreType())
    opts.SetScoreType(oefastrocs.OEShapeDatabaseType_Combo)
    print(opts.GetScoreType())

    ifs_0 = oechem.oemolistream()
    ifs_0.open(ref)

    ## Get size of reference set
    ref_size = 0
    for mol in ifs_0.GetOEMols():
        ref_size += 1

    ## Open reference file a second time
    ifs = oechem.oemolistream()
    ifs.open(ref)

    dbase = oefastrocs.OEShapeDatabase()
    moldb = oechem.OEMolDatabase()
    moldb.Open(ifs)

    dots = oechem.OEThreadedDots(10000, 200, 'conformers')
    dbase.Open(moldb, dots)

    qfs_0 = oechem.oemolistream()
    qfs_0.open(test)

    ## Get size of test set
    test_size = 0
    for mol in qfs_0.GetOEMols():
        test_size += 1

    ## Open test file a second time
    qfs = oechem.oemolistream()
    qfs.open(test)

    mcmol = oechem.OEMol()
    qmolidx = 0

    ## We use list of lists instead of numpy array because we want to store
    ## conformer information as a string
    output_matrix = np.zeros((ref_size, test_size))
    #output_matrix = [[0 for cols in range(query_size)] for rows in range(ref_size)]
    #conf_matrix = [[0 for cols in range(query_size)] for rows in range(ref_size)]

    ## Loop through all the compounds in the query set
    while oechem.OEReadMolecule(qfs, mcmol):
        ## Looping over each conf of the query molecule
        for q_index, q_conf in enumerate(mcmol.GetConfs()):
            #print('Reading conf {} of {}'.format(index+1, mcmol.NumConfs()))
            for score in dbase.GetSortedScores(q_conf, opts):
                dbmol_idx = score.GetMolIdx()
                #dbmol_confidx = score.GetConfIdx()
    	    #can change to just color or shape Tanimoto
                new_score = score.GetShapeTanimoto()
		#new_score = score.GetTanimotoCombo()
                if new_score > output_matrix[dbmol_idx][qmolidx]:
                    output_matrix[dbmol_idx][qmolidx] = new_score
                    #conf_matrix[qmolidx][dbmol_idx] = f'[{q_index}, {dbmol_confidx}]'
        qmolidx += 1

    ## Clean output matrix before saving information
    output_matrix = modify_3D_sim_mat(output_matrix, ref_groups, test_groups)   
    np.save(output, output_matrix)

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(description='''
    Calculate 3D similarity matrix of reference and test set
    ''', allow_abbrev=False)

    my_parser.add_argument('--ref',
            action='store',
            type=str,
            help='input .oeb database of reference SMILES',
            required=True)

    my_parser.add_argument('--test',
            action='store',
            type=str,
            help='input .oeb database of test SMILES',
            required=True)

    my_parser.add_argument('--ref_group',
            action='store',
            type=str,
            help='path to dictionary of flagged reference compound indices',
            required=True)

    my_parser.add_argument('--test_group',
            action='store',
            type=str,
            help='path to dictionary of flagged test compound indices',
            required=True)

    my_parser.add_argument('--output',
            action='store',
            type=str,
            help='name of output numpy array of similarity scores',
            required=True)

    args = my_parser.parse_args()
    calc_3D_sim(args.ref, args.test, args.ref_group, args.test_group, args.output)
