import sys
import os
import enhydra_modules
import multiprocessing
project_config = sys.argv[2] #path to project config file
code_config = sys.argv[1] #path to code config file

#try to open project config file, checks if it exists and can be read
try:
    fh_project_config = open(project_config, "r")
except OSError:
    print("Couldn't open the project configuration file, you may have to check if the name or permission of the file are correct")
    sys.exit()
#try to open code config file, checks if it exists and can be read
try:
    fh_code_config = open(code_config, "r")
except OSError:
    print("Couldn't open the config file, you may have to check if the name or permission of the file are correct")
    sys.exit()

parameters = enhydra_modules.read_config_file(fh_project_config, fh_code_config) #read the config files
fh_code_config.close()
fh_project_config.close()

enhydra_modules.check_parameters(parameters, code_config) #check if the parameters are correct

inputfiles = os.listdir(parameters['inputdir']) #read the input files

#create output directory
if os.path.isdir(parameters['outdir']): 
    exit("Outdir '%s' already exists, please change the parameter 'outdir'," % parameters['outdir'])
else:
    os.mkdir(parameters['outdir'])

#create log file
log_path = parameters['outdir'] + "/outlog"
outlog = open(log_path, "a")


print("Welcome to Enhydra\n")
outlog.write("Welcome to Enhydra\n")

#enhydra_modules.filter_length(parameters, inputfiles)
length_out_seq = parameters['outdir'] + "/length_filter"
length_out_stats = parameters['outdir'] + "/length_stats"
if not os.path.isdir(length_out_seq):
    os.mkdir(length_out_seq)
if not os.path.isdir(length_out_stats):
    os.mkdir(length_out_stats)

print ("Group filter step: length\n")
with multiprocessing.Pool(processes=parameters['max_process']) as pool:
      arg = [(parameters, outlog, files) for files in inputfiles]
      pool.starmap(enhydra_modules.filter_length, arg)

print ("Group filter step: min_species and anchor\n")

enhydra_modules.filter_groups(parameters, outlog)

print ("Alignment step\n")
enhydra_modules.run_mafft(parameters)

print ("Get alignment identities\n")
enhydra_modules.run_trimal(parameters)

print ("Making tables\n")
enhydra_modules.make_tables(parameters)