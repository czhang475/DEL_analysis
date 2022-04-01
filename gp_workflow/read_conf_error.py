import argparse
import numpy as np
import pickle

def read_error(logfile, indfile):
    '''
    Verify that conformers were properly generated for all molecules needing enumerated stereochemistry
    '''
    ## Parse error file to get indices of trouble molecules
    problem_ind = []
    with open(logfile) as fp:
        for index, line in enumerate(fp):
            if (index - 1) % 3 == 0:
                problem_ind.append(line.split(' ')[2])

    problem_ind = set([int(x) for x in problem_ind])

    ## Load in grouping dict and verify all trouble molecules had stereoisomers enumerated
    grouping = pickle.load(open(indfile, 'rb'))
    grouping_ind = set(grouping.keys())

    if problem_ind == grouping_ind:
        print('Conformers successfully generated for all molecules in {}'.format(logfile[5:-4]))
        return True
    else:
        print('Check molecules {}'.format(problem_ind ^ grouping_ind))
        return False



if __name__ == "__main__":
    my_parser = argparse.ArgumentParser(
            description="Returns indices of compounds that failed to generate conformers",
            allow_abbrev=False)

    my_parser.add_argument('--logfile',
            action='store',
            type=str,
            help='name of error file',
            required=True)

    my_parser.add_argument('--indfile',
            action='store',
            type=str,
            help='name of grouping dictionary',
            required=True)

    args = my_parser.parse_args()

    read_error(args.logfile, args.indfile)
