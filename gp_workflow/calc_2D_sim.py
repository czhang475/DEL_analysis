import numpy as np
import pandas as pd
import rdkit
import argparse
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

def gen_ecfp6(SMILES):
    set_mols = [Chem.MolFromSmiles(smi) for smi in SMILES]
    set_fps = [AllChem.GetMorganFingerprint(mol, 3, useCounts=True) for mol in set_mols]
    return set_fps

def ecfp6_tanimoto_matrix(row_fps, column_fps):
    matrix = np.zeros((len(row_fps), len(column_fps)))
    for i in range(len(row_fps)):
        for j in range(len(column_fps)):
            if DataStructs.TanimotoSimilarity(row_fps[i], column_fps[j]) > 0.1:
                matrix[i][j] = DataStructs.TanimotoSimilarity(row_fps[i], column_fps[j])
    return matrix

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(description='''
    Generate conformers given an input of SMILES strings
    ''', allow_abbrev=False)

    my_parser.add_argument('--ref',
            action='store',
            type=str,
            help='input .csv file of reference SMILES',
            required=True)

    my_parser.add_argument('--test',
            action='store',
            type=str,
            help='input .csv file of test SMILES',
            required=True)

    my_parser.add_argument('--output',
            action='store',
            type=str,
            help='name of output numpy array of similarity scores',
            required=True)

    args = my_parser.parse_args()

    ref_table = pd.read_csv(args.ref)
    test_table = pd.read_csv(args.test)

    ref_fps = gen_ecfp6(ref_table['SMILES'])
    test_fps = gen_ecfp6(test_table['SMILES'])

    output_matrix = ecfp6_tanimoto_matrix(ref_fps, test_fps)
    np.save(args.output, output_matrix)
