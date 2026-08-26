#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
import numpy as np
from os import path
from Bio import SeqIO
import os
import sys
from datetime import datetime
import shutil


DT_STRING = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def parse_args():
    parser = argparse.ArgumentParser(description='Extract sequences from reduced database')
    parser.add_argument('-f','--fasta', help='Fasta file containing the original database', required=True)
    parser.add_argument('-d','--db', help='TSV file containing the reduced database meta data', required=True)
    parser.add_argument('-p', '--prefix', help='Prefix. Default is input database')
    parser.add_argument('-o', '--out-dir', help='Output directory. Default is current working directory')
    parser.add_argument('-z', '--threads', help='Number of CPU threads to use. (default: %(default)s)', default='20')

    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    args = parser.parse_args()
    return args

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
      
def checkArgs(args):
    if args.prefix is None:
        prefix = path.basename(args.db)
        args.prefix = path.splitext(prefix)[0]
    if args.out_dir is None:
        file_dir= os.getcwd()
    else:
        file_dir = path.abspath(args.out_dir)
    args.out_dir = f'{file_dir}/cluster-extract-db_{args.prefix}_{DT_STRING}'
    if not path.exists(args.fasta):
        sys.exit(f'{args.fasta} not found')
    args.fasta = path.abspath(args.fasta)
    if not path.exists(args.db):
        sys.exit(f'{args.db} not found')
    args.db = path.abspath(args.db)
    args.threads = validateThreads(args.threads)
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    global RUN_LOG
    RUN_LOG = f'{scripts_path}/logs/cluster-extract-db_{DT_STRING}.log'
    if not path.exists(path.dirname(RUN_LOG)):
        os.mkdir(path.dirname(RUN_LOG))
    return args

def createOutputDir(args):
    dir_list = [args.out_dir]
    for item in dir_list:
        os.makedirs(item)
    printInfo(f'\nOutput directory: {args.out_dir}')
    os.chdir(args.out_dir)

def extractDB(args):
    db_df = pd.read_csv(args.db, sep='\t', header=0)
    random_set = set(db_df['accession'].to_list())
    records = []
    with open(args.fasta) as fopen:
        for rec in SeqIO.parse(fopen, 'fasta'):
            acc = rec.id
            if acc in random_set:
                records.append(rec)
    records_out = f'{args.out_dir}/{args.prefix}.fasta'
    file_out = open(records_out, 'w')
    SeqIO.write(records, file_out , 'fasta')
    printInfo(f'Created fasta: {records_out}')
    print(f'   ...{len(records)} records...')

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
    output_config = args.prefix+'_cluster-extract-db.config'
    
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

def main():
    args = parse_args()
    args = checkArgs(args)
    createConfig(args, var_dict={}, print_config=True) 
    createOutputDir(args)
    extractDB(args)
    output_config = createConfig(args, var_dict={})
    print('Done')

if __name__ == "__main__":
   main()       