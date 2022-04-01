from openeye import oechem, oeomega
import numpy as np
import pandas as pd
import tools
import argparse
import pickle

def gen_conf(infile, outfile, warnfile, flagfile):
    '''
    Generates conformers of the provided SMILES and outputs all molecules in an oeb file.
    '''
    ## Read in input .csv file containing all SMILES strings
    data = pd.read_csv(infile)

    ## Set up error catching
    errfs = oechem.oeofstream(warnfile)
    oechem.OEThrow.SetOutputStream(errfs)

    ## Initialize structures to store information
    new_table = pd.DataFrame(columns=['index', 'Molecule'])
    indices = []
    molecules = []

    ## Generate conformers for all SMILES strings
    for index, smi in enumerate(data['SMILES']):
        mol = tools.smiles_to_oemol(smi)
        oechem.OETriposAtomNames(mol)
        mol = tools.normalize_molecule(mol)
        cmol = None
        try:
            cmol = tools.generate_conformers(mol, max_confs=200)
            indices.append(index)
            molecules.append(cmol)

        ## if no specified stereochemistry, enumerate stereoisomers and generate conformers for all
        except Exception as e:
            oechem.OEThrow.Warning('Molecule {} returned an error\n{}'.format(index, str(e)))

        if cmol is None:
            for nmol in oeomega.OEFlipper(mol):
                oechem.OETriposAtomNames(nmol)
                try:
                    tools.generate_conformers(nmol, max_confs=200)
                    indices.append(index)
                    molecules.append(nmol)
                except Exception:
                    pass

    new_table['index'] = indices
    new_table['Molecule'] = molecules

    ## Save entire dataframe of generated conformers into database
    tools.write_dataframe_to_file(new_table, outfile)

    ## Record the compounds with enumerated stereoisomers and their indices for later analysis
    groupings = {}
    group = []
    for index, value in enumerate(indices[:-1], 1):
        # start sorting indices if the value of entry is same as the previous
        if indices[index] == indices[index-1]:
            group.append(index-1)
            # edge case to close off group if it occurs at the very end of the list
            if index == len(indices)-1:
                group.append(index)
                groupings[value] = group
        # stop storing indices once value has changed, close off list, append it and make a new list
        elif indices[index] != indices[index-1] and group:
            group.append(index-1)
            groupings[value] = group
            group = []
        else:
            continue

    pickle.dump(groupings, open(flagfile, 'wb'))

    ## Read error -- will output in std.out (slurm-{jobnumber}.out) if there are issues with any conformers
    tools.read_error(warnfile, flagfile)

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(
            description="Generate conformers given an input of SMILES strings",
            allow_abbrev=False)

    my_parser.add_argument('--infile',
            action='store',
            type=str,
            help='input .csv file of SMILES',
            required=True)

    my_parser.add_argument('--outfile',
            action='store',
            type=str,
            help='name of output .oeb of generated conformers',
            required=True)


    my_parser.add_argument('--warnfile',
            action='store',
            type=str,
            help='file of warnings for compounds that failed to generate conformers',
            required=True)

    my_parser.add_argument('--flagfile',
            action='store',
            type=str,
            help='dictionary containing compound indices and their stereoisomers',
            required=True)

    args = my_parser.parse_args()

    gen_conf(args.infile, args.outfile, args.warnfile, args.flagfile)
