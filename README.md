"### this is the introduction of the usage of cmCluster code" 
"### here we strat cmCluster at the count matrix from cellranger output files"

#cd the file folder where all of these scripts are
$cd script_file_folder

#run create_cmCluster_sh.py to get the shell file to do cmCluster
$python3 create_cmCluster_sh.py -i input_file_folder -o output_file_folder -w 7

#run shell file(run_cmCluster.sh) create in your output_file_folder to finish cmCluster
$sh run_cmCluster.sh

#the result of cmCluster will be put in output_file_folder
