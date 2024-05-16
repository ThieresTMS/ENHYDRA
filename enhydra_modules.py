import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import numpy as np
import multiprocessing

def read_config_file(fh_project, fh_code):
    parameters = {}
    parameters['inputdir'] = read_line("inputdir", fh_project)
    parameters['outdir'] = read_line("outdir", fh_project)
    parameters['min_species'] = int(read_line("min_species", fh_project))
    parameters['anchor'] = read_line("anchor", fh_project)
    parameters['mafft'] = read_line("mafft", fh_code)
    parameters['trimal'] = read_line("trimal", fh_code)
    parameters['max_process'] = int(read_line("max_process", fh_project))
    return parameters


def read_line(parameter, file):
    for line in file:
        line = line.rstrip()
        match = re.search("^%s" % parameter, line)
        if match:
            line = re.sub('^\s*%s\s*=\s*' % parameter, '', line)
            line = re.sub('#.*$', '', line)
            line = re.sub('\s*$', '', line)
            return line


def check_parameters(parameters, code_config):
    if parameters['inputdir'] == "":
        exit("No path to the inputdir was specified in your project's configuration file, please fill the parameter 'inputdir'.")
    if not os.path.isdir(parameters['inputdir']):
        exit("The path to your project's species files isn't a valid file, please check if the path in 'infile' is correct: %s" % parameters['inputdir'])
    if not os.access(parameters['inputdir'], os.R_OK):
        exit("You don't have permission to read in your project's species file: %s, please redefine your permissions." % parameters['inputdir'])

    if parameters['outdir'] == "":
        exit("No path to outdir was specified in project_config, please open this file and fill the parameter 'outdir'")
    if parameters['min_species'] == "":
        exit("No minimum number of species per group was specified in your project's configuration file, please fill the the parameter 'min_species'.")

    if parameters['anchor'] == "":
        exit("No anchor species was specified in project_config, please open this file and fill the parameter 'anchor'")
    if parameters['max_process'] == "":
        exit("Maximum number of process was not specified in project_config, please open this file and fill the parameter 'max_process'")

    if parameters['trimal'] == "":
        exit("No path to trimal was specified in code_config at %s, please open this file and fill the parameter 'trimal'" % code_config)
    if not os.path.isfile(parameters['trimal']):
        exit("The executable of trimal wasn't found in the specified path, please check if the path is correct: %s" % parameters['trimal'])
    if not os.access(parameters['trimal'], os.R_OK):
        exit("You don't have permission to execute the trimal file specified at code_config, please check permissions or replace the file")

    if parameters['mafft'] == "":
        exit("No path to mafft was specified in code_config at %s, please open this file and fill the parameter 'mafft'" % code_config)
    if not os.path.isfile(parameters['mafft']):
        exit("The executable of mafft wasn't found in the specified path, please check if the path is correct: %s" % parameters['mafft'])
    if not os.access(parameters['mafft'], os.R_OK):
        exit("You don't have permission to execute the mafft file specified at code_config, please check permissions or replace the file")

def filter_length(parameters, outlog, file):
    inputfile = parameters['inputdir'] + "/" + file
    if os.stat(inputfile).st_size == 0:
        return
    outfile_s_path = parameters['outdir'] + "/length_stats" + "/" + file + "_lengthstats"
    outfile_f_path = parameters['outdir'] + "/length_filter" + "/" + file + "_lengthfilter"
    outstats = open(outfile_s_path, "a")
    outfile = open(outfile_f_path, "a")
    lengths = []
    length_sum = 0
    length_data = {}
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        seq = Seq(seq_record.seq)
        seq_id = seq_record.id
        lenght = len(seq)
        length_data[seq_id] = lenght
        lengths.append(lenght)
        length_sum = length_sum + lenght
    total_seqs = lengths
    mean = statistics.mean(lengths)
    median = statistics.median(lengths)
    stddev = 0
    if len(lengths) > 1:
        stddev = statistics.stdev(lengths)
    outstats.write("##Overal sequence length stats\n")
    outstats.write("Total seqs: %s\nAverage: %s\nMedian: %s\nSD: %s\n" % (total_seqs, mean, median, stddev))
    outstats.write("##Sequence lengths (sorted form smallest to largest)\n")
    outstats.write("#SequenceID\tLength\tPercentageDifFromAvg\n")
    keys = list(length_data.keys())
    values = list(length_data.values())
    sorted_value_index = np.argsort(values)
    length_data_sorted = {keys[i]: values[i] for i in sorted_value_index}
    for key in length_data_sorted:
        dist = length_data_sorted[key]/mean
        outstats.write("%s\t%s\t%s\n" % (key, length_data_sorted[key], dist))
    outstats.close()
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        seq = Seq(seq_record.seq)
        greater = (mean + (2*stddev))
        smaller = (mean - (2*stddev))
        if (len(seq) < smaller) or (len(seq) > greater):
            outlog.write("Sequence %s in group %s removed for length filter step\n" % (seq_record.id, file))
        else:
            outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def filter_groups(parameters, outlog):
    out = parameters['outdir'] + "/group_filter"
    if not os.path.isdir(out):
        os.mkdir(out)
    in_dir_path = parameters['outdir'] + "/length_filter"
    files = os.listdir(in_dir_path)
    for file in files:
        file_fields = file.split(".")
        path_to_file = in_dir_path + "/" + file
        outfile_path = out + "/" + file
        species_ids = []
        for seq_record in SeqIO.parse(path_to_file, "fasta"):
            ids_fields = seq_record.id.split("|")
            sid = ids_fields[0]
            species_ids.append(sid)
        uniq_ids = list(set(species_ids))
        if parameters['anchor'] not in uniq_ids:
            outlog.write("Group %s not contain any sequence of anchor species %s. Group removed from analysis\n" % (file_fields[0], parameters['anchor']))
            continue
        if len(uniq_ids) < int(parameters['min_species']):
            outlog.write("Number of species in group %s is less than minimum required %s , group removed from analysis\n" % (file_fields[0], parameters['min_species']))
        else:
            os.system('cp %s %s' % (path_to_file, outfile_path))


def run_mafft(parameters):

    dirpath = parameters['outdir'] + '/group_filter'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        output_seq = outdir + "/" + file + ".aln"
        os.system('%s --auto --quiet --thread %s %s > %s' % (parameters['mafft'], parameters['max_process'], input_seq, output_seq))


def run_trimal(parameters):
    dirpath = parameters['outdir'] + '/alignment'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/ident_alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        outident = outdir + "/" + file + ".ident"
        os.system('%s -sident -in %s > %s' % (parameters['trimal'], input_seq, outident)),


def make_tables(parameters):
    tables = parameters['outdir'] + "/tables"
    if not os.path.isdir(tables):
        os.mkdir(tables)
    ortho_mean = {}
    dirpath1 = parameters['outdir'] + "/alignment"
    input_files1 = os.listdir(dirpath1)
    group2anchor_path = tables + "/group2anchor.tsv"
    group2anchor= open(group2anchor_path, "a")
    dirpath2 = parameters['outdir'] + "/ident_alignment"
    input_files2 = os.listdir(dirpath2)
    group2mean_path = tables + "/group2mean.tsv"
    anchor2mean_path = tables + "/anchor2mean.tsv"
    anchor2mean = open(anchor2mean_path, "a")
    group2mean = open(group2mean_path, "a")
    for file2 in input_files2:
        aux2 = file2.split(".")
        group_name2 = aux2[0]
        ident_files = dirpath2 + '/' + file2
        ident_file = open(ident_files, "r")
        for line in ident_file:
            line = line.rstrip()
            match = re.search("identity:", line)
            if match:
                values = line.split()
                mean_percent = values[-1]
                ortho_mean[group_name2] = mean_percent
                group2mean.write("%s\t%s\n" % (group_name2, mean_percent))
    for file1 in input_files1:
        aux = file1.split(".")
        group_name = aux[0]
        seq_file = dirpath1 + '/' + file1
        for seq_record in SeqIO.parse(seq_file, "fasta"):
            ids_fields = seq_record.id.split("|")
            specie = ids_fields[0]
            sequence_id = ids_fields[1]
            if specie == parameters['anchor']:
                anchor2mean.write("%s\t%s\n" % (sequence_id, ortho_mean[group_name]))
                group2anchor.write("%s\t%s\n" % (group_name, sequence_id))
                    