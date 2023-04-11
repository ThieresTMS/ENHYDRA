import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import numpy as np


def read_config_file(fh_project, fh_code):
    parameters = {}
    parameters['inputdir'] = read_line("inputdir", fh_project)
    parameters['outdir'] = read_line("outdir", fh_project)
    parameters['min_species'] = read_line("min_species", fh_project)
    parameters['mafft'] = read_line("mafft", fh_code)
    parameters['trimmal'] = read_line("trimmal", fh_code)
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


##inserir aqui a checagem para o mafft##

##inserir aqui a checagem para o trimmal##


def filter_length(parameters, inputfiles):
    length_out_seq = parameters['outdir'] + "/length_filter"
    length_out_stats = parameters['outdir'] + "/length_stats"
    if not os.path.isdir(length_out_seq):
        os.mkdir(length_out_seq)
    if not os.path.isdir(length_out_stats):
        os.mkdir(length_out_stats)
    for file in inputfiles:
        inputfile = parameters['inputdir'] + "/" + file
        if os.stat(inputfile).st_size == 0:
            continue
        outfile_s_path = length_out_stats + "/" + file + "_lengthstats"
        outfile_f_path = length_out_seq + "/" + file + "_lengthfilter"
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
        outstats.write("##Sequence lenghts (sorted form smallest to largest)\n")
        outstats.write("#SequenceID\tLenght\tPercentageDifFromAvg\n")
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
                continue
                #print("%s é menor que %s ou é maior que %s" % (len(seq), smaller, greater))
            else:
                outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def filter_groups(parameters, outlog):
    out = parameters['outdir'] + "/group_filter"
    if not os.path.isdir(out):
        os.mkdir(out)
    in_dir_path = parameters['outdir'] + "/length_filter"
    files = os.listdir(in_dir_path)
    for file in files:
        path_to_file = in_dir_path + "/" + file
        outfile_path = out + "/" + file
        species_ids = []
        for seq_record in SeqIO.parse(path_to_file, "fasta"):
            ids_fields = seq_record.id.split("|")
            sid = ids_fields[0]
            species_ids.append(sid)
        uniq_ids = list(set(species_ids))
        if len(uniq_ids) < int(parameters['min_species']):
            outlog.write("Number of species in group %s is less than minimum required %s , group removed from analysis\n" % (file, parameters['min_species']))
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
        os.system('%s --auto  %s > %s' % (parameters['mafft'], input_seq, output_seq ) )


def run_trimmal(parameters):
    dirpath = parameters['outdir'] + '/alignment'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/ident_alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        outident = outdir + "/" + file + ".ident"
        os.system('%s -sident -in %s > %s' % (parameters['trimmal'], input_seq, outident))

