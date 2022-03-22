#! /usr/bin/bash
#SBATCH --job-name=pipeline
#SBATCH --partition=titanx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --gres=gpu:titan:1

#Function to call to run the actual code
source activate oepython

script_path=/DFS-L/DATA/mobley/chrisyz1/anagenex/Paper/scripts
cd $script_path

# Initialize variables
reftype='reference_bbs'
testtype='test_bbs'
seed=12

# Create new directory to store all files
dirname=data/split_$seed
mkdir -p $dirname

# Create reference and test split and save files in new directory
input_files=$dirname/input_files
mkdir -p $input_files
python reference_test_split.py --seed ${seed} --path "${input_files}"

# Generate conformers for reference and test sets
conf_files=$dirname/conformers
mkdir -p $conf_files

warning_files=$dirname/warnings
mkdir -p $warning_files

python gen_conf.py --infile "${input_files}/${reftype}_${seed}.csv" --outfile "${conf_files}/${reftype}_${seed}.oeb" --warnfile "${warning_files}/${reftype}_${seed}.log" --flagfile "${warning_files}/${reftype}_${seed}.pkl"
python gen_conf.py --infile "${input_files}/${testtype}_${seed}.csv" --outfile "${conf_files}/${testtype}_${seed}.oeb" --warnfile "${warning_files}/${testtype}_${seed}.log" --flagfile "${warning_files}/${testtype}_${seed}.pkl"

# Confirm that conformers are generated for all reference and test inputs
#python read_conf_error.py --logfile "${dirname}/${reftype}_${seed}.log" --indfile "${dirname}/${reftype}_${seed}.pkl"
#python read_conf_error.py --logfile "${dirname}/${testtype}_${seed}.log" --indfile "${dirname}/${testtype}_${seed}.pkl"

# Generate 2D and 3D similarity matrices
mat_files=$dirname/similarity_matrices
mkdir -p $mat_files

python calc_2D_sim.py --ref "${input_files}/${reftype}_${seed}.csv" --test "${input_files}/${testtype}_${seed}.csv" --output "${mat_files}/${reftype}_${testtype}_${seed}_2D.npy"
python calc_3D_sim.py --ref "${conf_files}/${reftype}_${seed}.oeb" --test "${conf_files}/${testtype}_${seed}.oeb" --ref_group "${warning_files}/${reftype}_${seed}.pkl" --test_group "${warning_files}/${testtype}_${seed}.pkl" --output "${mat_files}/${reftype}_${testtype}_${seed}_3D.npy"

# Remove extra compound stereoisomers from 3D similarity matrix
#python clean_3D_sim_matrix.py --matrix "${mat_files}/${reftype}_${testtype}_${seed}_3D.npy" --ref_group "${warning_files}/${reftype}_${seed}.pkl" --test_group "${warning_files}/${testtype}_${seed}.pkl"

rank_files=$dirname/rankings
mkdir -p $rank_files
# Return an ordered list of which test compounds are predicted to be active
python get_rank.py --matrix "${mat_files}/${reftype}_${testtype}_${seed}_2D.npy" --test "${input_files}/${testtype}_${seed}.csv" --output "${rank_files}/${reftype}_${testtype}_${seed}_2D_ranked.csv"
python get_rank.py --matrix "${mat_files}/${reftype}_${testtype}_${seed}_3D.npy" --test "${input_files}/${testtype}_${seed}.csv" --output "${rank_files}/${reftype}_${testtype}_${seed}_3D_ranked.csv" 
