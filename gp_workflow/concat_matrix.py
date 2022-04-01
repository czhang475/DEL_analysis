import numpy as np
import pandas as pd
from tools import *
from os import listdir
from os.path import isfile, join
import argparse

def concat_columns(refpath, testpath):
    # Get size of reference conformer set
    ref_size = len(read_file_to_dataframe(refpath))
    col_size = len([f for f in listdir(testpath) if isfile(join(testpath, f))])
    full_mat = np.zeros((ref_size, col_size))

    # Read in each column and save to full_mat
    for i in range(col_size):
        column = np.load(testpath+'/test_set_{}_3D.npy'.format(i))
        column_flat = [item for sublist in column for item in sublist]
        full_mat[:,i] = column_flat

    np.save(testpath+'/test_mat_3D.npy', full_mat)

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(description='''Concatenate columns of a similarity matrix''')
    my_parser.add_argument('--refpath', action='store', type=str, help='path to reference conformer file')
    my_parser.add_argument('--testpath', action='store', type=str, help='path to matrix column directory')

    args = my_parser.parse_args()
    concat_columns(args.refpath, args.testpath)
