import sys
import os
import enhydra_modules
import multiprocessing
project_config = sys.argv[2]
code_config = sys.argv[1]

try:
    fh_project_config = open(project_config, "r")
except OSError:
    print("Couldn't open the project configuration file, you may have to check if the name or permission of the file are correct")
    sys.exit()
try:
    fh_code_config = open(code_config, "r")
except OSError:
    print("Couldn't open the config file, you may have to check if the name or permission of the file are correct")
    sys.exit()

parameters = code_modules.read_config_file(fh_project_config, fh_code_config)
fh_code_config.close()
fh_project_config.close()

code_modules.check_parameters(parameters, code_config)

inputfiles = os.listdir(parameters['inputdir'])

if os.path.isdir(parameters['outdir']):
    exit("Outdir '%s' already exists, please change the parameter 'outdir'," % parameters['outdir'])
else:
    os.mkdir(parameters['outdir'])

log_path = parameters['outdir'] + "/outlog"
outlog = open(log_path, "a")
print("Welcome to Code")
outlog.write("Welcome to Code\n")

#code_modules.filter_length(parameters, inputfiles)
length_out_seq = parameters['outdir'] + "/length_filter"
length_out_stats = parameters['outdir'] + "/length_stats"
if not os.path.isdir(length_out_seq):
    os.mkdir(length_out_seq)
if not os.path.isdir(length_out_stats):
    os.mkdir(length_out_stats)

with multiprocessing.Pool(processes=parameters['max_process']) as pool:
      arg = [(parameters, files) for files in inputfiles]
      pool.starmap(code_modules.filter_length, arg)

code_modules.filter_groups(parameters, outlog)

code_modules.run_mafft(parameters)

code_modules.run_trimal(parameters)

code_modules.make_tables(parameters)
