#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
from os import path
import os
import sys
from datetime import datetime
from natsort import natsort_keygen


DT_STRING = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

import functions_primer_check as pc

def parse_args():
    parser = argparse.ArgumentParser(description='Primer check using BLASTn')
    parser.add_argument('-f','--file', help='TSV file containing primers. Must have columns "[primer_name, seq, amp_name, amp_len, oligo_type (left, right, probe)]"', required=True)
    parser.add_argument('-p', '--prefix', help='Output prefix. Default uses tsv file name')
    parser.add_argument('-x', '--off-target-db', help='Path to off-target blast database. If provided, primers will be checked for potential off-targets', default=None, type=str)
    parser.add_argument('-y', '--target-db', help='Path to target blast database. If provided, primers will be checked whether they can potentially catch all targets', default=None, type=str)
    parser.add_argument('--target-mismatch', help='Max number of mismatches when comparing to target blast db (default: %(default)s).', default=3, type=int)
    # parser.add_argument('--target-size-diff', help='+/- bp for on-target size (default: %(default)s).', default=30, type=int)
    parser.add_argument('-o', '--out-dir', help='Output directory. Default is <current directory>/primer-check_<prefix>_<date>')
    parser.add_argument('-z', '--threads', help='Number of CPU threads to use. (default: %(default)s)', default='10')
    parser.add_argument('-c', '--on-target-column', help='Prefix for on-target column names. (default: %(default)s)', default='on_target')
    parser.add_argument('-g', '--off-target-column', help='Prefix for off-target column names. (default: %(default)s)', default='off_target')
    parser.add_argument('--allowance-3-prime', help="Number of mismatches allowed in the first 3 bases of the 3' end for off-target (default: %(default)s)", default=0, type=int)


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

def checkArgs(args):
    if args.prefix is None:
        prefix = path.basename(args.file)
        args.prefix = prefix.replace(f'.{prefix.split('.')[-1]}', '')
    if args.out_dir is None:
        file_dir= os.getcwd()
    else:
        file_dir = path.abspath(args.out_dir)
    args.out_dir = f'{file_dir}/primer-check_{args.prefix}_{DT_STRING}'
    if not path.exists(args.file):
        sys.exit(f'Primer tsv file {args.file} not found')
    args.file = path.abspath(args.file)
    args.threads = validateThreads(args.threads)
    args.off_target_db = pc.retrieveBlastDB(args.off_target_db)
    args.target_db = pc.retrieveBlastDB(args.target_db)
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    global RUN_LOG
    RUN_LOG = f'{scripts_path}/logs/primer-check_{DT_STRING}.log'
    if not path.exists(path.dirname(RUN_LOG)):
        os.mkdir(path.dirname(RUN_LOG))
    return args

def createOutputDir(args):
    dir_list = [args.out_dir]
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
    output_config = args.prefix+'_primer-check.config'
    
    if print_config:
        printInfo(tab)
    else:
        printInfo(f'\nCreated {output_config}')
        tab.to_csv(output_config, sep='\t', index=False, header=True)
    return output_config

################# Primers #################
def readPrimers(args):
    df1 = pd.read_csv(args.file, sep='\t', header=0)
    columns = ['primer_name', 'seq','amp_name', 'amp_len', 'oligo_type']
    for col in columns:
        if col not in df1.columns:
            sys.exit(f'{col} not found in provided primer_tsv\nProvide the following columns:\n[primer_name, seq, amp_name, amp_len, oligo_type (left, right, probe)]')
    return df1


################# Main #################


def main():
    args = parse_args()
    args = checkArgs(args)
    printMessage('Primer check using BLASTn')
    createConfig(args, var_dict={}, print_config=True) 
    createOutputDir(args)
    createConfig(args, var_dict={}) 
    primer_df = readPrimers(args)
    print(primer_df)
    blast_query  = f'{args.out_dir}/{args.prefix}_query.fasta'
    pc.createBlastQuery(primer_df.copy(), blast_query)
 
    if args.off_target_db is not None:
        n_seqs = pc.targetCount(args.off_target_db)
        blast_out = f'{args.out_dir}/{args.prefix}_blast_off-targets.tsv'
        offtarget_dict = pc.searchOffTargets(primer_df, n_seqs, args.allowance_3_prime, blast_query, args.threads, args.off_target_db, blast_out)
        primer_df[f'{args.off_target_column}_no'] = primer_df['primer_name'].map(offtarget_dict['off_targets_no']).astype(str) + '/' + str(n_seqs)
        primer_df[args.off_target_column] = primer_df['primer_name'].map(offtarget_dict['off_targets'])
        primer_df[f'{args.off_target_column}_species'] = primer_df['primer_name'].map(offtarget_dict['off_target_species'])
        print(primer_df)
        

    if args.target_db is not None:
        n_seqs = pc.targetCount(args.target_db)
        blast_out = f'{args.out_dir}/{args.prefix}_blast_on-targets.tsv'
        target_dict = pc.searchTargets(primer_df, n_seqs, blast_query, args.threads, args.target_db, blast_out, args.target_mismatch)
        primer_df[args.on_target_column] = primer_df['primer_name'].map(target_dict['success'])
        primer_df[args.on_target_column] = primer_df[args.on_target_column].astype(str) + '/' + str(n_seqs)
        primer_df[f'{args.on_target_column}_fail_species'] = primer_df['primer_name'].map(target_dict['fail_species']) 
        primer_df[f'{args.on_target_column}_fails'] = primer_df['primer_name'].map(target_dict['fails']) 
    
    print(primer_df)
    primer_out = f'{args.out_dir}/{args.prefix}_primer-check.tsv'
    primer_df.to_csv(primer_out, sep='\t', index=False)

    print('Done')

if __name__ == "__main__":
   main()   


