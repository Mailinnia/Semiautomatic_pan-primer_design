#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
import numpy as np
from os import path
import os
import sys
from datetime import datetime
import shutil
import random
from Bio import AlignIO
from Bio.Seq import Seq
import math
import itertools
import multiprocessing
from natsort import natsort_keygen
import string
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product

# imported varvamp scripts
import varvamp_alignment
import functions_primer_check as pc


DT_STRING = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def parse_args():
    parser = argparse.ArgumentParser(description='Primer design using varVamp')
    parser.add_argument('-a','--aln', help='Alignment file', required=True)
    parser.add_argument('-t','--thresholds', help='Consensus thresholds (default: "%(default)s")', default='0.75 0.80 0.85', type=str)
    parser.add_argument('-n','--max-ambigs', help='List of max number of ambiguous characters in a primer  (default: "%(default)s")', default='3 4 5 6', type=str)
    parser.add_argument('-m','--mode', help='varVamp mode (default: %(default)s). Automatically switches to single mode if no primers are found for qpcr.', default='qpcr', choices=['qpcr', 'single'])
    parser.add_argument('-p', '--prefix', help='Output prefix. Default uses alignment name')
    parser.add_argument('-b', '--primer-prefix', help='Primer prefix (default: %(default)s).', default='Primer', type=str)
    parser.add_argument('-x', '--off-target-db', help='Path to off-target blast database. If provided, primers will be checked for potential off-targets', default=None, type=str)
    parser.add_argument('--allowance-3-prime', help="Number of mismatches allowed in the first 3 bases of the 3' end for off-target (default: %(default)s)", default=0, type=int)
    parser.add_argument('-y', '--target-db', help='Path to target blast database. If provided, primers will be checked whether they can potentially catch all targets', default=None, type=str)
    parser.add_argument('--target-mismatch', help='Max number of mismatches when comparing to target blast db (default: %(default)s).', default=3, type=int)
    parser.add_argument('--opt-amp', help='Optimal amplicon length when single mode is used (default: %(default)s).', default='150', type=str)
    parser.add_argument('--max-amp', help='Max amplicon length when single mode is used (default: %(default)s).', default='300', type=str)
    parser.add_argument('--min-len', help='Min primer length (default: %(default)s).', default=18, type=int)
    parser.add_argument('--opt-len', help='Optimal primer length (default: %(default)s).', default=21, type=int)
    parser.add_argument('--max-len', help='Max primer length (default: %(default)s).', default=24, type=str)
    parser.add_argument('--min-qamp', help='Min q-amplicon length when qPCR mode is used (default: %(default)s).', default=70, type=int)
    parser.add_argument('--max-qamp', help='Max q-amplicon length when qPCR mode is used (default: %(default)s).', default=200, type=int)
    parser.add_argument('--min-temp', help='Min primer temp (default: %(default)s).', default=56, type=int)
    parser.add_argument('--max-temp', help='Max primer temp (default: %(default)s).', default=63, type=int)
    parser.add_argument('--opt-temp', help='Optimal primer temp (default: %(default)s).', default=60, type=int)
    parser.add_argument('-o', '--out-dir', help='Output directory. Default is <current directory>/primer-design_<prefix>_<date>')
    parser.add_argument('-z', '--threads', help='Number of CPU threads to use. (default: %(default)s)', default='10')
    parser.add_argument('-r', '--report-title', help='Title of report. Default is prefix')

    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    args = parser.parse_args()
    return args

def runProcess(cmd, pipe=None, ignore_returncode=False, shell=False):
    base_names = []
    for item in cmd:
        if item.endswith('/'):
            item_split = item.split('/')
            basename = f'{item_split[len(item_split)-2]}/'
        else:
            basename = path.basename(item)
        base_names.append(basename)
    # base_names= [path.basename(item) for item in cmd]
    print_cmd = f"\n\t{' '.join(base_names)}"
    print(print_cmd)
    log_out = suppressLog()
    log_out.write(f'\n{" ".join(cmd)}\n')
    log_out.close()
    if shell:
        print('\tRunning in shell')
        cmd = ' '.join(cmd)
        print(cmd)
        process = subprocess.run(cmd, stderr = suppressLog(), shell=True)
    elif pipe is None:
        process = subprocess.run(cmd, stderr = suppressLog())
    else:
        process = subprocess.run(cmd, stdout=pipe)
    if ignore_returncode: 
        return process.returncode
    elif process.returncode != 0:
        fail_msg = f'{cmd} FAILED'
        quitLog(fail_msg)

################# Logs and messages #################
def printMessage(message):
    white_len = (90-4-len(message))/2
    if not white_len.is_integer():
        white_len += 0.5
    if not white_len %2 == 0:
        white_len +=1 
    white_len = int(white_len)
    white_space = ' '*white_len
    pound_len = white_len*2+4+len(message)
    before = f'\n{"#"*pound_len}\n'
    mid = f'##{white_space*2}{" "*len(message)}##\n'
    after = f'{"#"*pound_len}\n'
    final_message = f'{before}{mid}##{white_space}{message}{white_space}##\n{mid}{after}'
    print(final_message)
    log_out = suppressLog()
    log_out.write(f'\n{final_message}\n')
    log_out.close()

def printInfo(message):
    msg = f'\n[INFO] {message}'
    print(msg)
    log_out = suppressLog()
    log_out.write(f'\n{msg}\n')
    log_out.close()

def printWarning(message):
    msg = f'\n[WARNING] {message}'
    print(msg)
    log_out = suppressLog()
    log_out.write(f'\n{msg}\n')
    log_out.close()

def quitLog(fail_msg):
    fail_log = RUN_LOG
    if path.isfile(fail_log):
        fout = open(fail_log, 'a')
    else:
        fout = open(fail_log, 'w')
    print(f'\n{fail_msg}')
    fout.write(fail_msg)
    fout.close()
    print('\nQuitting pipeline')
    quit()

def suppressLog():
    suppress_log = RUN_LOG
    if path.isfile(suppress_log):
        fout = open(suppress_log, 'a')
    else:
        fout = open(suppress_log, 'w')
    return fout

def validateThreads(threads):
    cpu_no = os.cpu_count()
    if int(threads) > cpu_no:
        printInfo('There are only', cpu_no, 'threads available.')
        threads = str(cpu_no)
    return threads

def retrieveBlastDB(db):
    if db is not None:
        if path.isdir(db):
            found = False
            for file in os.listdir(db):
                if file.endswith('.ndb'):
                    db_name = path.basename(file).replace('.ndb', '')
                    found=True
                    break
            if found:
                db = f'{path.abspath(db)}/{db_name}'
            else:
                sys.exit(f'Cannot find blast db files in {db}')
        else:
            sys.exit(f'{db} is not a directory')
    return db

def targetCount(args):
    db_path = path.dirname(args.target_db)
    for file in os.listdir(db_path):
        if file.endswith('.njs'):
            njs = path.basename(file)
            njs_path = f'{db_path}/{njs}'
            break

    with open(njs_path, 'r') as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('"number-of-sequences":'):
                n_seqs = int(line.split(':')[1].strip().strip(','))
                return n_seqs

def parseThresholds(thresholds):
    thresholds_split = thresholds.split()
    thresh_list = []
    for thresh in thresholds_split:
        thresh_float = float(thresh)
        if thresh_float == 0 or thresh_float >1:
            sys.exit(f'{thresh} is not a valid value. Thresholds must be >0 and <=1')
        else:
            thresh_list.append(str(thresh_float))
    return sorted(thresh_list)

def checkArgs(args):
    if args.prefix is None:
        prefix = path.basename(args.aln)
        args.prefix = prefix.replace(f'.{prefix.split(".")[-1]}', '')
    if args.report_title is None:
        args.report_title = args.prefix
    if args.out_dir is None:
        file_dir= os.getcwd()
    else:
        file_dir = path.abspath(args.out_dir)
    args.out_dir = f'{file_dir}/primer-design_{args.prefix}_{DT_STRING}'
    if not path.exists(args.aln):
        sys.exit(f'Alignment {args.aln} not found')
    args.aln = path.abspath(args.aln)
    args.threads = validateThreads(args.threads)
    args.thresholds = parseThresholds(args.thresholds)
    args.max_ambigs = sorted([x.strip() for x in args.max_ambigs.split()])[::-1] #reverse list
    args.off_target_db = retrieveBlastDB(args.off_target_db)
    args.target_db = retrieveBlastDB(args.target_db)
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    args.rmd = f'{scripts_path}/primer_design.Rmd'
    args.rscript = f'{scripts_path}/primer_design.R'
    args.varvamp_dir = f'{args.out_dir}/varvamp'
    args.post_dir = f'{args.out_dir}/post_processing'
    global RUN_LOG
    RUN_LOG = f'{scripts_path}/logs/primer-design_{DT_STRING}.log'
    if not path.exists(path.dirname(RUN_LOG)):
        os.mkdir(path.dirname(RUN_LOG))
    global VARVAMP_CONFIG
    VARVAMP_CONFIG =f'{scripts_path}/varvamp_config.py'
    return args

def createOutputDir(args):
    dir_list = [args.out_dir, args.varvamp_dir, args.post_dir]
    for item in dir_list:
        os.makedirs(item)
    printInfo(f'\nOutput directory: {args.out_dir}')
    os.chdir(args.out_dir)

################# Config #################
def createConfig(args, var_dict, print_config=False):  
    args_tuple = []
    for arg in vars(args):
        argument = arg
        value = getattr(args, arg) 
        if value is None:
            value = 'None'
        args_tuple.append((argument, value))
   
    for var in var_dict.keys():
        args_tuple.append((var, var_dict[var]))
    args_tuple.append(('run_log', RUN_LOG))
    args_tuple.append(('date', DT_STRING))
    tab = pd.DataFrame(args_tuple)
    tab.columns = ['Argument', 'Value']
    output_config = args.prefix+'_primer-design.config'
    
    if print_config:
        printInfo(tab)
    else:
        printInfo(f'\nCreated {output_config}')
        tab.to_csv(output_config, sep='\t', index=False, header=True)
    return output_config

def add2varDict(key_list, var_list, var_dict=None):
    if var_dict is None:
        var_dict = {}
    for key, var in zip(key_list, var_list):
        var_dict[key]=var
    return var_dict

def executeR(args, output_config):
    printMessage('GENERATING HTML REPORT')
    cmd = ['Rscript', args.rscript]
    R_cmd = ['-c', output_config]
    full_cmd = cmd+R_cmd
    runProcess(full_cmd)

################# varVamp #################
def parseVarVampLog(outdir):
    log = []
    with open(f'{outdir}/varvamp_log.txt') as fp:
        for line in fp:
            line = line.strip()
            if len(line) > 1:
                log.append(line)
    return log[-1]

def modifyConfig(args):
    '''
    Modifies varvamp default config with new parameters
    '''
    config_out = f'{args.out_dir}/varvamp_config.py'
    modify = False
    if args.min_len != 18 or args.opt_len != 21 or args.max_len != 24 or args.max_qamp != 200 or args.min_qamp != 70 or args.min_temp != 56 or args.max_temp != 63 or args.opt_temp != 60:
        modify =True
        fout = open(config_out, 'w')
        with open(VARVAMP_CONFIG, 'r') as fp:
            for line in fp:
                if line.startswith('PRIMER_SIZES'):
                    line = line.replace('(18, 24, 21)', f'({args.min_len}, {args.max_len}, {args.opt_len})')
                if line.startswith('QAMPLICON_LENGTH'):
                    line = line.replace('(70, 200)', f'({args.min_qamp}, {args.max_qamp})')
                if line.startswith('PRIMER_TMP'):
                    line = line.replace('(56, 63, 60)', f'({args.min_temp}, {args.max_temp}, {args.opt_temp})')
                if line.startswith('QPROBE_TMP'):
                    line = line.replace('(64, 70, 67)', f'({args.min_temp+7}, {args.max_temp+7}, {args.opt_temp+7})')
                fout.write(line)
        fout.close()
        
    return modify, config_out

def callVarVamp(args, threshold, ambig, mode):
    outdir = f'{args.varvamp_dir}/{mode}-t{threshold}-a{ambig}'
    name = f'{args.prefix}-t{threshold}-a{ambig}'
    if mode == 'qpcr':
        cmd = ['varvamp', 'qpcr', '-th', args.threads, '--name', name, '-a', ambig, '-t', threshold, args.aln, outdir]
    else:
        cmd = ['varvamp', 'single', '-th', args.threads, '--name', name, '-a', ambig, '-t', threshold, '-ol', args.opt_amp, '-ml', args.max_amp, args.aln, outdir]
    
    modify, config_file = modifyConfig(args)
    if modify:
        printInfo('Modifying varvamp config file')
        varvamp_config_call = [f'VARVAMP_CONFIG={config_file}']
        cmd = varvamp_config_call+cmd
        return_code = runProcess(cmd, ignore_returncode=True, shell=True)
    else:
        return_code = runProcess(cmd, ignore_returncode=True)
    if return_code !=0:
        primers_found =False
        printInfo(f'varVamp failed for threshold {threshold} with max ambiguity {ambig} in {mode} mode: {parseVarVampLog(outdir)}')           
    else:
        primers_found = True
    return outdir, primers_found

def avoidInsanity(mode, tuple_list, min_threshold, max_ambigs):
    """
    Avoid using parameters we know will fail. If primers are not found at lower
    thresholds, they will not be found at higher ones. If ambig value does not work
    at lower thresholds, it will for sure not work at higher ones. Therefore we 
    adjust the ambig values and check if it makes sense to continue in qPCR mode, and
    whether it makes sense to continue increasing threshold.
    """
    df1 = pd.DataFrame(tuple_list)
    df1.columns = ['mode', 'result', 'threshold', 'max_ambig', 'dir']
    df2 = df1[(df1['result'] == 'pass') & (df1['threshold'] == min_threshold)].copy()
    skip = False
    if df2.empty:
        printWarning(f'No primers found for threshold {min_threshold}. It does not make sense to try higher thresholds. Will skip.')
        skip = True # doesn't make sense to try higher thresholds
    modes = df2['mode'].to_list()
    if mode == 'qpcr':
        if 'qpcr' in modes:
            mode = 'qpcr'
        else:
            printWarning(f'qpcr not successful for {min_threshold}. It does not make sense to try qpcr mode for higher thresholds. Will switch permanently to single mode')
            mode = 'single'
    else:
        mode = mode
    ambigs = sorted(list(set(df2['max_ambig'].to_list())))[::-1]
    if ambigs != max_ambigs:
        if len(ambigs) != 0:
            printWarning(f'Adjusting max_ambigs from {max_ambigs} to {ambigs}')
            
    return mode, ambigs, skip

def printSkipAmbigs(max_ambigs, ambig):
    remaning_ambigs = max_ambigs[max_ambigs.index(ambig)+1:]
    if len(remaning_ambigs) > 0:
        printWarning(f'Skipping remaining ambigs {remaning_ambigs}')

def findPrimers(args):
    """
    Calls varVamp for each threshold and ambig combination. If qpcr mode fails,
    then it tries single mode
    """
    printInfo('Finding primers with varVamp')
    outdirs = []
    min_threshold = args.thresholds[0]
    thresh_mode = args.mode
    max_ambigs = args.max_ambigs
    for threshold in args.thresholds:
        #sanity check
        if threshold != min_threshold:
            thresh_mode, max_ambigs, skip = avoidInsanity(thresh_mode, outdirs, min_threshold, max_ambigs)
            min_threshold = threshold
            if skip:
                #add each threshold we've yet to loop through
                index = args.thresholds.index(threshold) # find current place in list
                for i in range(index, len(args.thresholds)):
                    outdirs.append((thresh_mode, 'skip', args.thresholds[i], pd.NA, pd.NA))
                break
        ambig_mode = thresh_mode
        for ambig in max_ambigs:
            outdir, primers_found = callVarVamp(args, threshold, ambig, ambig_mode)
            if primers_found:
                outdirs.append((ambig_mode, 'pass', threshold, ambig, outdir))
            else:
                outdirs.append((ambig_mode, 'fail', threshold, ambig, pd.NA))
                shutil.rmtree(outdir)
                if ambig_mode == 'qpcr':
                    printInfo('Trying single mode')
                    outdir, primers_found = callVarVamp(args, threshold, ambig, 'single')
                    if primers_found:
                        outdirs.append(('single', 'pass', threshold, ambig, outdir))
                        printWarning(f'Switching to single mode for {threshold}')
                        ambig_mode = 'single' #if qpc failed for current ambig, it makes no sense to try at lower ambigs.
                    else:
                        shutil.rmtree(outdir)
                        outdirs.append(('single', 'fail', threshold, ambig, pd.NA))
                        printSkipAmbigs(max_ambigs, ambig)
                        break #if current ambig fails, then skip lower ones as they will definitely fail
                else:
                    printSkipAmbigs(max_ambigs, ambig)
                    break #if current ambig fails, then skip lower ones as they will definitely fail
    df1 = pd.DataFrame(outdirs)
    df1.columns = ['mode', 'result', 'threshold', 'max_ambig', 'dir']
    df1['threshold'] = df1['threshold'].astype(float)
    df1['max_ambig'] = df1['max_ambig'].astype('Int64')
    df1_out = f'{args.post_dir}/{args.prefix}_status.tsv'
    df1.to_csv(df1_out, sep='\t', index=False)
    tab = df1.drop(columns=['dir']).copy()
    print(tab)
    return df1, df1_out


################# Check primers #################
def checkPrimers(args, primer_df):
    blast_query  = f'{args.post_dir}/{args.prefix}_query.fasta'
    pc.createBlastQuery(primer_df.copy(), blast_query)
 
    if args.off_target_db is not None:
        n_seqs = pc.targetCount(args.off_target_db)
        blast_out = f'{args.post_dir}/{args.prefix}_blast_off-targets.tsv'
        offtarget_dict = pc.searchOffTargets(primer_df, n_seqs, args.allowance_3_prime, blast_query, args.threads, args.off_target_db, blast_out)
        primer_df['off_target_species'] = primer_df['primer_name'].map(offtarget_dict['off_target_species'])
        primer_df['off_targets'] = primer_df['primer_name'].map(offtarget_dict['off_targets'])
    else:
        primer_df['off_target_species']  = 'n.d.'
        primer_df['off_targets'] = 'n.d.'

    if args.target_db is not None:
        n_seqs = pc.targetCount(args.target_db)
        blast_out = f'{args.post_dir}/{args.prefix}_blast_off-targets.tsv'
        target_dict = pc.searchTargets(primer_df, n_seqs, blast_query, args.threads, args.target_db, blast_out, args.target_mismatch)
        primer_df['on_targets'] = primer_df['primer_name'].map(target_dict['success'])
        primer_df['on_targets'] = primer_df['on_targets'].astype(str) + '/' + str(n_seqs)
        primer_df['on_target_fail_species'] = primer_df['primer_name'].map(target_dict['fail_species'])
        primer_df['on_target_fails'] = primer_df['primer_name'].map(target_dict['fails']) 
    else:
        primer_df['on_targets']  = 'n.d.'
        primer_df['on_target_fail_species'] = 'n.d.'
        primer_df['on_target_fails'] = 'n.d.'

    return primer_df
    

################# Process varVamp outputs #################
def potentialRegions(args, varvamp_df):
    """
    Merges all the different region bed files to one dataframe for use in plot
    """ 
    printInfo('Collating all potential primer regions')
    pass_df = varvamp_df[varvamp_df['result']=='pass']
    pass_thresholds = sorted(list(set(pass_df['threshold'].to_list())))
    df2_list = []
    for threshold in pass_thresholds:
        threshold_df = pass_df[pass_df['threshold']==threshold].copy()
        varvamp_dirs = threshold_df['dir'].to_list()
        modes = threshold_df['mode'].to_list()
        df1_list = []   
        for mode, varvamp_dir in zip(modes, varvamp_dirs):
            ##Parsing potential primer regions file
            df1 = pd.read_csv(f'{varvamp_dir}/data/primer_regions.bed', sep='\t', header=None)
            df1.columns = ['chr', 'start', 'end', 'name']
            df1 = df1.drop(['chr', 'name'],axis=1)
            df1['data_type']='primer_region'
            df1_list.append(df1)
            if mode == 'qpcr':
                df1 = pd.read_csv(f'{varvamp_dir}/data/probe_regions.bed', sep='\t', header=None)
                df1.columns = ['chr', 'start', 'end', 'name']
                df1 = df1.drop(['chr', 'name'],axis=1)
                df1['data_type']='probe'
                df1_list.append(df1)
            else:
                df1 = pd.read_csv(f'{varvamp_dir}/data/all_primers.bed', sep='\t', header=None)
                df1.columns = ['chr', 'start', 'end', 'name', 'penalty', 'strand']
                df1['data_type']=df1['name'].apply(lambda x: x.split('_')[0])
                df1 = df1.drop(['chr', 'name', 'penalty', 'strand'],axis=1)
                df1_list.append(df1)

        df2 = pd.concat(df1_list, ignore_index=True).drop_duplicates()
        df2['threshold']=threshold
        df2_list.append(df2)
    regions_df = pd.concat(df2_list)
    regions_out = f'{args.post_dir}/{args.prefix}_potential_regions.tsv'
    regions_df.to_csv(regions_out, sep='\t', index=False)
    print(regions_df)
    # return regions_out
    return regions_df

def collatePrimers(varvamp_df):
    """
    Merges all passed primers to one dataframe
    """
    printInfo('Collating all primers')
    pass_df = varvamp_df[varvamp_df['result']=='pass']
    pass_thresholds = sorted(list(set(pass_df['threshold'].to_list())))[::-1] #decreasing order
    primers_list = []
    for threshold in pass_thresholds:
        threshold_df = pass_df[pass_df['threshold']==threshold].copy()
        modes = threshold_df['mode'].to_list()
        max_ambigs = threshold_df['max_ambig'].to_list()
        varvamp_dirs = threshold_df['dir'].to_list()
        primer_list = []   
        for mode, max_ambig, varvamp_dir in zip(modes, max_ambigs, varvamp_dirs):
            ##Parsing potential primer regions file
            if mode == 'single':
                df1 = pd.read_csv(f'{varvamp_dir}/primers.tsv', sep='\t', header=0)
                df1.columns = ['amp_name', 'amp_len', 'primer_name', 'primer_name_all_primers', 'pool', 'start', 'stop', 'seq', 'size', 'gc_best', 'temp_best', 'mean_gc', 'mean_temp', 'penalty', 'off_target_amplicons']
                df1['oligo_type'] = df1['primer_name'].copy().apply(lambda x: x.split('_')[-1])
                primers = df1[['amp_name', 'amp_len', 'primer_name', 'seq', 'size', 'start', 'stop', 'gc_best', 'temp_best', 'mean_gc', 'mean_temp', 'penalty', 'off_target_amplicons', 'oligo_type']].copy()
            else:
                df1 = pd.read_csv(f'{varvamp_dir}/qpcr_primers.tsv', sep='\t', header=0)
                df1.columns = ['qpcr_scheme', 'oligo_type', 'start', 'stop', 'seq', 'size', 'gc_best', 'temp_best', 'mean_gc', 'mean_temp', 'penalty', 'off_target_amplicons']
                df1 = df1.rename(columns={'qpcr_scheme':'amp_name'})
                df1['primer_name']=df1['amp_name']+'_'+df1['oligo_type']
                amplicon_lengths = df1[df1['oligo_type'] == 'LEFT'].merge(df1[df1['oligo_type'] == 'RIGHT'], on='amp_name', suffixes=('_left', '_right'))
                amplicon_lengths['amp_len'] = amplicon_lengths['stop_right'] - amplicon_lengths['start_left']
                amplicon_lengths = amplicon_lengths[['amp_name', 'amp_len']]
                df2 = df1.merge(amplicon_lengths, on='amp_name')
                primers = df2[['amp_name', 'amp_len', 'primer_name', 'seq', 'size', 'start', 'stop', 'gc_best', 'temp_best', 'mean_gc', 'mean_temp', 'penalty', 'off_target_amplicons', 'oligo_type']].copy()
            primers['max_ambig'] = int(max_ambig)
            primers['mode'] = mode
            primer_list.append(primers)
        primer_df = pd.concat(primer_list, ignore_index=True).drop_duplicates()
        primer_df['threshold']= float(threshold)
        primers_list.append(primer_df)
    primers_df = pd.concat(primers_list,ignore_index=True).reset_index(drop=True)
    return primers_df

def filterAmplicons(prefix,primers_df):
    """
    Checks whether primersets are duplicated. Assigns same name to duplicated primer sets.
    """
    df1 = primers_df[primers_df['oligo_type'].isin(['LEFT','RIGHT'])]
    df2 = df1[['threshold', 'max_ambig', 'amp_name', 'oligo_type', 'seq']]
    df3 = df2.pivot(index=['threshold', 'max_ambig', 'amp_name'], columns='oligo_type', values='seq').reset_index()
    natsort_key = natsort_keygen()
    df3 = df3.sort_values(by=['threshold', 'max_ambig', 'amp_name'], ascending=[False, True, True], key=lambda col: col.map(natsort_key) if col.name == 'amp_name' else col).reset_index(drop=True)
    df3.columns.name = None 
    df3['new_amp_name'] = pd.NA
    first_occurrences = {}
    i = 0
    for idx, row in df3.iterrows():
        #create key tuples {(left, right):amp_name}
        key = (row['LEFT'], row['RIGHT'])
        if key not in first_occurrences:
            i += 1
            if len(str(i)) ==1:
                primer_val = f'0{i}'
            else:
                primer_val = i
            new_name = f'{prefix}_{primer_val}'
            first_occurrences[key] =  new_name
            df3.at[idx, 'new_amp_name'] = new_name
        else:
            #check if tuple key is in dict. If in dict, match the name to the first occurence
            new_name = first_occurrences[key]
            df3.at[idx, 'new_amp_name'] = new_name
    naming_dict = df3.set_index('amp_name')['new_amp_name'].to_dict()
    filter_primers = df1.copy()
    filter_primers ['new_amp_name'] = filter_primers['amp_name'].map(naming_dict)
    filter_primers['probe_suffix'] = pd.NA
    return filter_primers, naming_dict

def incrementLetter(label):
    """
    Increment a label like 'a' -> 'b', ..., 'z' -> 'aa'.
    """
    if not label:
        return 'a'
    if label[-1] != 'z':
        return label[:-1] + chr(ord(label[-1]) + 1)
    return incrementLetter(label[:-1]) + 'a'

def filterProbes(df1, naming_dict):
    """
    Checks whether probes are duplicates. If the case should happen
    that a primerset has more than one probe (due to ambig)
    then probe is assigned a, b, etc. 
    """
    df2 = df1[df1['oligo_type']=='PROBE'].copy()
    df2['new_amp_name'] = df2['amp_name'].map(naming_dict)
    df3 = df2.copy()
    df3['probe_suffix']= pd.NA
    match_dict = {}
    for idx, row in df3.iterrows():
        #create key tuples {(left, right):amp_name}
        key = row['new_amp_name']
        seq = row['seq']
        if key not in match_dict:
            match_dict[key] =  {'seq_set':{seq}, 'suffix':'a', 'first_status':'not_changed', 'first_idx':idx}
        else:
            #check if tuple key is in dict. If in dict, match the name to the first occurence
            cur_suffix = match_dict[key]['suffix']
            if seq not in match_dict[key]['seq_set']:
                match_dict[key]['seq_set'].add(seq)
                new_suffix = incrementLetter(cur_suffix)
                df3.at[idx, 'probe_suffix'] = new_suffix
                match_dict[key]['suffix'] = new_suffix
                if match_dict[key]['first_status']=='not_changed':
                    df3.at[match_dict[key]['first_idx'], 'probe_suffix'] = cur_suffix
                    match_dict[key]['first_status'] = 'changed' 
    return df3

def filterPrimers(args, varvamp_df):
    """
    Filter primer and probes for duplicates and reassign names.
    Output new primer table
    """
  
    primers_df = collatePrimers(varvamp_df)
    print('\t... filtering duplicated primer sets ...')
    amplicon_filter, naming_dict = filterAmplicons(args.primer_prefix,primers_df.copy())
    probe_filter = filterProbes(primers_df.copy(), naming_dict)

    ## merge filtered dataframe 
    df1 = pd.concat([amplicon_filter, probe_filter], ignore_index=True)
    natsort_key = natsort_keygen()
    df1 = df1.sort_values(by=['threshold', 'max_ambig', 'amp_name', 'oligo_type'], ascending=[False, True, True, True], key=lambda col: col.map(natsort_key) if col.name == 'amp_name' else col).reset_index(drop=True)

    ## rename primers to F, R, Pr for LEFT, RIGHT, PROBE
    df1['new_type']=df1['oligo_type'].map({'LEFT':'F', 'RIGHT':'R', 'PROBE':'Pr'}) 
    df1['new_primer_name']=df1['new_amp_name']+df1['probe_suffix'].fillna('')+'_'+df1['new_type']
    ambig = ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B']
    df1['ambig_count'] = df1['seq'].apply(lambda seq: sum(seq.count(char) for char in ambig)) # count how many ambiguities in seq

    ## Need dictionary to filter amplicon regions
    filter_dict = df1.set_index('primer_name')[['new_primer_name', 'new_amp_name']].to_dict() 

    ## select columns for output and rename
    df2 = df1[['new_primer_name', 'seq', 'new_amp_name', 'amp_len', 'oligo_type', 'size', 'start', 'stop', 'ambig_count', 'gc_best', 'temp_best', 'mean_gc', 'mean_temp', 'penalty', 'threshold', 'max_ambig', 'mode']].copy()
    df2['oligo_type']=df2['oligo_type'].str.lower()
    df2 = df2.rename(columns={'new_primer_name':'primer_name', 'new_amp_name':'amp_name', 'size':'primer_size'})

    df2['note'] = df2.copy().groupby('primer_name')['threshold'].transform(lambda x: f"Also in threshold(s) {', '.join(map(str, sorted(x.unique())))}" if len(x.unique()) > 1 else pd.NA)
    df2['note'] = df2.apply(lambda row: row['note'].replace(f", {row['threshold']}", "").replace(f"{row['threshold']}, ", "").replace(f"{row['threshold']}", "").strip(", ") if pd.notna(row['note']) else row['note'], axis=1)
    df2 = df2.sort_values(by=['threshold', 'max_ambig', 'amp_name', 'oligo_type'], ascending=[False, True, True, True], key=lambda col: col.map(natsort_key) if col.name == 'amp_name' else col).reset_index(drop=True)
    df2_out = f'{args.post_dir}/{args.prefix}_primers_unfilt.tsv'
    df2.to_csv(df2_out, sep='\t', index=False) # output unfiltered collated primers
    filter_primers = df2.drop_duplicates(subset=['primer_name'], keep='first').reset_index(drop=True)
    filter_primer = checkPrimers(args, filter_primers)
    filter_out = f'{args.out_dir}/{args.prefix}_primers.tsv'
    filter_primers.to_csv(filter_out, sep='\t', index=False) # output filtered primers - this is the main file output
    print(filter_primers)
    return filter_out, filter_dict

def ampliconRegions(args, varvamp_df, filter_dict):
    """
    Merges the primerset regions into one dataframe.
    Then filters the primersets for those that have been 
    kept after checking for duplicates
    """ 
    printInfo('Collating amplicon regions')
    pass_df = varvamp_df[varvamp_df['result']=='pass']
    pass_thresholds = sorted(list(set(pass_df['threshold'].to_list())))
    primer_regions_list = []
    for threshold in pass_thresholds:
        threshold_df = pass_df[pass_df['threshold']==threshold].copy()
        varvamp_dirs = threshold_df['dir'].to_list()
        # max_ambigs = threshold_df['max_ambig'].to_list()
        primer_region_list = []   
        for varvamp_dir in varvamp_dirs:
        # for varvamp_dir, ambig in zip(varvamp_dirs, max_ambigs):
            ##Parsing primer regions file
            file = f'{varvamp_dir}/primers.bed'
            df1 = pd.read_csv(file, sep='\t', header=None, comment="#" )
            df1.columns = ['chr', 'start', 'end', 'name', 'pool', 'strand', 'seq']
            df2 = df1.copy().drop(['chr', 'pool', 'seq'], axis=1)
            df2['new_name'] = df2['name'].map(filter_dict['new_primer_name'])
            df2['amp_name'] = df2['name'].map(filter_dict['new_amp_name'])
            # df2['max_ambig'] = ambig
            primer_region_list.append(df2)
        primer_region_df = pd.concat(primer_region_list, ignore_index=True).drop_duplicates()
        primer_region_df['threshold']=threshold
        primer_regions_list.append(primer_region_df)

    df3 = pd.concat(primer_regions_list, ignore_index=True)

    # adding oligo type for easier filtering downstream
    df3['oligo_type']=df3['name'].copy().apply(lambda x: x.split('_')[-1].lower()) 

    # removing records which are not in filter_dict
    primer_region = df3[df3['new_name'].notna()].reset_index(drop=True)

    # Retrieving left start and right end to get regions for full amplicon
    df4 = primer_region.copy()
    df4 = df4[df4['oligo_type'].isin(['left', 'right'])].reset_index(drop=True) #don't need probes
    df4 = df4.drop_duplicates(subset=['amp_name', 'threshold', 'oligo_type'], keep='first').reset_index(drop=True)

    # Pivot without dropping duplicate amp_name entries
    df4 = (df4.pivot(index=['amp_name', 'threshold'], columns='oligo_type', values=['start', 'end']).reset_index())

    # Flatten the multi-level columns
    df4.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col for col in df4.columns]

    #create amplicon regions
    df5 = df4[['start_left','end_right','amp_name','threshold']].copy()
    df5 = df5.rename(columns={'start_left':'start', 'end_right':'end', 'amp_name':'name'})
    df5['oligo_type']='amplicon'

    #select primer set regions
    df6 = primer_region[['start', 'end', 'new_name', 'oligo_type', 'threshold']].copy()
    df6 = df6.rename(columns={'new_name':'name'})

    #merge regions
    region_df = pd.concat([df6, df5], ignore_index=True)
    natsort_key = natsort_keygen() #natural sort of names
    region_df = region_df.sort_values(by=['threshold', 'start', 'name'], key=lambda col: col.map(natsort_key) if col.name == 'name' else col).reset_index(drop=True)
    region_df = region_df.rename(columns={'oligo_type':'data_type'}).drop_duplicates(keep='first').reset_index(drop=True)
    region_out = f'{args.post_dir}/{args.prefix}_amplicon_regions.tsv'
    region_df.to_csv(region_out, sep='\t', index=False)
    # print(region_df)
    # return region_out
    return region_df

################# Recreate entropy data #################
def entropy(chars, states):
    """
    input is a list of characters or numbers.
    calculate relative shannon's entropy. relative values are
    achieved by using the number of possible states as the base
    """
    ent = 0
    n_chars = len(chars)
    # only one char is in the list
    if n_chars <= 1:
        return ent
    # calculate the number of unique chars and their counts
    value, counts = np.unique(chars, return_counts=True)
    probs = counts/n_chars
    if np.count_nonzero(probs) <= 1:
        return ent

    for prob in probs:
        ent -= prob*math.log(prob, states)

    return ent

def alignmentEntropy(args, threshold):
    """
    calculate the entropy for every position in an alignment.
    return pandas df.
    """
    print(f'\t... preprocessing alignment ...')
    preprocessed_alignment = varvamp_alignment.preprocess(args.aln)
    # clean alignment of gaps
    print(f'\t... cleaning alignment ...')
    # alignment_cleaned, gaps_to_mask = varvamp_alignment.process_alignment(preprocessed_alignment, float(threshold), int(args.threads))
    alignment_cleaned, gaps_to_mask = varvamp_alignment.process_alignment(preprocessed_alignment, float(threshold))
    ent_tuple = []
    # iterate over alignment positions and the sequences
    for nuc_pos in range(0, len(alignment_cleaned[0][1])):
        pos = []
        for seq_number in range(0, len(alignment_cleaned)):
            pos.append(alignment_cleaned[seq_number][1][nuc_pos])
        ent_tuple.append((nuc_pos, entropy(pos, 4)))
    # create df
    ent_df = pd.DataFrame(ent_tuple)
    ent_df.columns=['pos', 'ent']
    ent_df["mean"] = ent_df["ent"].rolling(10, center=True).mean()
    return ent_df

def entropyData(args):
    printInfo('Calculating entropies')
    print('\tThis may be slow if there are many unique gaps in the alignment')
    ent_list = []
    for threshold in args.thresholds:
        print(f'\t## {threshold} ##')
        df1 = alignmentEntropy(args, threshold)
        df1['threshold']=float(threshold)
        ent_list.append(df1)
    ent_df = pd.concat(ent_list, ignore_index=True)
    ent_out = f'{args.post_dir}/{args.prefix}_entropy.tsv'
    ent_df.to_csv(ent_out, sep='\t', index=False)
    print(ent_df)
    # return ent_out
    return ent_df

################# Create static images #################
def entropySubplot(prefix,ax, df1, thresh):
    """
    creates the entropy subplot
    """
    print('\t\t... adding entropies ...')
    df2 = df1[df1['threshold']==thresh]
    # - create entropy df
    ax[0].fill_between(df2["pos"], df2["ent"], color="gainsboro", label="Entropy")
    ax[0].plot(df2["pos"], df2["mean"], color="black", label="Rolling mean entropy", linewidth=0.5)
    ax[0].set_ylim((0, 1))
    ax[0].set_xlim(0, max(df2["pos"]))
    ax[0].set_ylabel("Normalized Shannon's entropy")
    ax[0].set_title(f"Consensus threshold {thresh}\n{prefix}")
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    #
    ## ax is mutable, and therefore doesn't need to be returned
    return df2['pos'].max()

def regionSubplot(ax, primer_df, thresh, location, color, description):
    """
    Creates the primer regions subplot using a DataFrame.
    """
    primer_df = primer_df[primer_df['threshold']==thresh].reset_index(drop=True)
    #
    if not primer_df.empty:
        # Plot each region as a horizontal line
        locations = [location] * len(primer_df) 
        ax[1].hlines(locations, primer_df['start'], primer_df['end'], linewidth=5, color=color)
        # Add a single legend entry (using the first row to avoid redundant labels)
        ax[1].hlines(location, primer_df.iloc[0]['end'], primer_df.iloc[0]['end'], 
                    label=description, linewidth=5, color=color)
    ## ax is mutable, and therefore doesn't need to be returned

def allRegionSubplot(ax, primer_df, thresh):
    """
    Creates the primer regions subplot using a DataFrame.
    """
    print('\t\t... adding primers ...')
    primer_df = primer_df[primer_df['threshold']==thresh]
    left_df = primer_df[primer_df['data_type']=='LEFT'].copy().reset_index(drop=True)
    right_df = primer_df[primer_df['data_type']=='RIGHT'].copy().reset_index(drop=True)
    left_color = 'dimgrey'
    right_color = 'darkgrey'
    loc = 0.85
    if not left_df.empty:
    # Plot each region as a horizontal line
        locations = [loc]*len(left_df)
        ax[1].hlines(locations, left_df['start'], left_df['end'], linewidth=5, color=left_color)
        # Add a single legend entry (using the first row to avoid redundant labels)
        ax[1].hlines(loc, left_df.iloc[0]['end'], left_df.iloc[0]['end'], 
                    label='All forward primers', linewidth=5, color=left_color)
        loc = 0.80
    if not right_df.empty:
        locations = [loc]*len(right_df)
        ax[1].hlines(locations, right_df['start'], right_df['end'], linewidth=5, color=right_color)
        # Add a single legend entry (using the first row to avoid redundant labels)
        ax[1].hlines(loc, right_df.iloc[0]['end'], right_df.iloc[0]['end'], 
                    label='All reverse primers', linewidth=5, color=right_color)
        loc = 0.75
    ## ax is mutable, and therefore doesn't need to be returned
    return loc
   
def ampliconSubplot(ax, primer_df, thresh, loc, max_pos):
    """
    Creates the primer regions subplot using a DataFrame.
    """
    print('\t\t... adding amplicons ...')
    primer_df = primer_df[primer_df['threshold']==thresh]
    #
    amplicon_df = primer_df[primer_df['data_type']=='amplicon'].copy().sort_values(by=['start', 'end']).reset_index(drop=True)
    amplicon_df['text']=amplicon_df['name'].apply(lambda x: x.split('_')[-1])
    amplicon_df['text_loc'] = amplicon_df['end']-((amplicon_df['end']-amplicon_df['start'])/2)
    #
    probe_df = primer_df[primer_df['data_type']=='probe'].copy().reset_index(drop=True)
    #
    prim_df = primer_df[primer_df['data_type'].isin(['left', 'right'])].copy().reset_index(drop=True)
    #
    amplicon_color = '#1f77b4'
    probe_color = '#00C5CD'
    prim_color = '#d62728'
    # Plot each region as a horizontal line
    locations = [loc]*len(amplicon_df)
    ax[1].hlines(locations, amplicon_df['start'], amplicon_df['end'], linewidth=5, color=amplicon_color)
    # Add text annotations
    txt_location = loc-0.15  # Fixed y-coordinate for text
    ratio_val = round(max_pos/300, ndigits=0)
    prev_start = -ratio_val
    prev_end = -ratio_val
    min_loc = []
    for idx, row in amplicon_df.iterrows():
        if (prev_start - ratio_val <= row['start'] <= prev_start + ratio_val) or (prev_end - ratio_val <= row['end'] <= prev_end + ratio_val):
            # print('Imma gonna yshift')
            txt_location -= 0.06
            min_loc.append(txt_location)
        else:
            txt_location = loc-0.15
        prev_start = row['start']
        prev_end = row['end'] 
        ax[1].text(row['text_loc'], txt_location, row['text'], fontsize=8, ha='center', va='bottom')
    #
    # Add a single legend entry (using the first row to avoid redundant labels)
    ax[1].hlines(loc, amplicon_df.iloc[0]['end'], amplicon_df.iloc[0]['end'], 
                 label='Amplicons', linewidth=5, color=amplicon_color)   
    # Plot each region as a horizontal line
    locations = [loc]*len(prim_df)
    ax[1].hlines(locations, prim_df['start'], prim_df['end'], linewidth=5, color=prim_color)
    # Add a single legend entry (using the first row to avoid redundant labels)
    ax[1].hlines(loc, prim_df.iloc[0]['end'], prim_df.iloc[0]['end'], 
                 label='Primers', linewidth=5, color=prim_color)
    if not probe_df.empty:
        # Plot each region as a horizontal line
        locations = [loc-0.05]*len(probe_df)
        ax[1].hlines(locations, probe_df['start'], probe_df['end'], linewidth=5, color=probe_color)
        # Add a single legend entry (using the first row to avoid redundant labels)
        ax[1].hlines(loc-0.05, probe_df.iloc[0]['end'], probe_df.iloc[0]['end'], 
                    label='Probes', linewidth=5, color=probe_color)
    if len(min_loc) ==0:
        min_loc = loc
    else:
        min_loc = min(min_loc)
    return min_loc
    ## ax is mutable, and therefore doesn't need to be returned

def varVampPlot(args, entropy_df, regions_df, amplicon_df):
    printInfo('Creating plots')
    for threshold in args.thresholds:
        threshold = float(threshold)
        print(f'\t{threshold}')
        fig_out = f'{args.post_dir}/{args.prefix}_t{threshold}.svg'
        fig, ax = plt.subplots(2, 1, figsize=[22, 6], squeeze=True, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        fig.subplots_adjust(hspace=0)

        if not regions_df[regions_df['threshold']==threshold].empty:
            primer_regions = regions_df[regions_df['data_type']=='primer_region'].copy()
            probe_regions = regions_df[regions_df['data_type']=='probe'].copy()
            left_right_regions = regions_df[regions_df['data_type'].isin(['LEFT', 'RIGHT'])].copy()

            max_pos = entropySubplot(args.report_title, ax, entropy_df.copy(), threshold)
            print('\t\t... adding regions...')
            regionSubplot(ax, primer_regions, threshold, location=0.95, color="darkorange",description="Possible primer regions")
            print('\t\t... adding probes ...')
            regionSubplot(ax, probe_regions, threshold, location=0.90, color="#8B1C62",description="Possible probe regions" )
            loc = allRegionSubplot(ax, left_right_regions, threshold)
            min_loc = ampliconSubplot(ax, amplicon_df, threshold, loc, max_pos)
            # finalize
            print('\t\t... finalizing ...')
            ax[1].spines['right'].set_visible(False)
            ax[1].spines['left'].set_visible(False)
            ax[1].spines['bottom'].set_visible(False)
            ax[1].axes.get_yaxis().set_visible(False)
            ax[1].set_xlabel("Alignment position")
            ax[1].set_ylim((min_loc-0.15, 1))
            fig.legend(loc=(0.835, 0.6))

            # save fig
            print('\t\t... attempting to save ...')
            fig.savefig(fig_out, bbox_inches='tight')
            print('\t\t... successful! ...')
            plt.close()

################# Main #################

def main():
    args = parse_args()
    args = checkArgs(args)
    printMessage('Primer design using varVamp')
    createConfig(args, var_dict={}, print_config=True) 
    createOutputDir(args)
    createConfig(args, var_dict={}) 
    entropy_df = entropyData(args)
    varvamp_df, status_tsv = findPrimers(args)
    regions_df= potentialRegions(args, varvamp_df)
    primers_tsv, filter_dict = filterPrimers(args, varvamp_df)
    amplicon_df = ampliconRegions(args, varvamp_df, filter_dict)
    varVampPlot(args, entropy_df, regions_df, amplicon_df)
    output_config = createConfig(args, var_dict={})
    executeR(args, output_config)
    print('Done')

if __name__ == "__main__":
   main()   


# dfs = dict()
# for file in os.listdir('./'):
#     if file.endswith('entropy.tsv'):
#         dfs['entropy']=file
#     elif file.endswith('potential_regions.tsv'):
#         dfs['regions']=file
#     elif file.endswith('amplicon_regions.tsv'):
#         dfs['amplicon']=file

# regions_df = pd.read_csv(dfs['regions'], sep='\t', header=0)
# amplicon_df = pd.read_csv(dfs['amplicon'], sep='\t', header=0)
# entropy_df = pd.read_csv(dfs['entropy'], sep='\t', header=0)
# args_post_dir = os.getcwd()
# args_prefix = dfs['entropy'].replace('_entropy.tsv', '')
# args_report_title = 'Orthoflavivirus'
# threshold = '0.75'

# for threshold in ['0.75','0.80','0.85']:
#     threshold = float(threshold)
#     if not regions_df[regions_df['threshold']==threshold].empty:
#         print(f'\t{threshold}')
#         fig_out = f'{args_post_dir}/{args_prefix}_t{threshold}.svg'
#         fig, ax = plt.subplots(2, 1, figsize=[22, 6], squeeze=True, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
#         fig.subplots_adjust(hspace=0)
#         #
#         primer_regions = regions_df[regions_df['data_type']=='primer_region'].copy()
#         probe_regions = regions_df[regions_df['data_type']=='probe'].copy()
#         left_right_regions = regions_df[regions_df['data_type'].isin(['LEFT', 'RIGHT'])].copy()
#         #
#         max_pos = entropySubplot(args_report_title, ax, entropy_df.copy(), threshold)
#         print('\t\t... adding regions...')
#         regionSubplot(ax, primer_regions, threshold, location=0.95, color="darkorange",description="Possible primer regions")
#         print('\t\t... adding probes ...')
#         regionSubplot(ax, probe_regions, threshold, location=0.90, color="#8B1C62",description="Possible probe regions" )
#         loc = allRegionSubplot(ax, left_right_regions, threshold)
#         min_loc = ampliconSubplot(ax, amplicon_df, threshold, loc, max_pos)
#         # finalize
#         print('\t\t... finalizing ...')
#         ax[1].spines['right'].set_visible(False)
#         ax[1].spines['left'].set_visible(False)
#         ax[1].spines['bottom'].set_visible(False)
#         ax[1].axes.get_yaxis().set_visible(False)
#         ax[1].set_xlabel("Alignment position")
#         ax[1].set_ylim((min_loc-0.15, 1))
#         fig.legend(loc=(0.835, 0.6))
#         #
#         # save fig
#         print('\t\t... attempting to save ...')
#         fig.savefig(fig_out, bbox_inches='tight')
#         print('\t\t... successful! ...')
#         plt.close()