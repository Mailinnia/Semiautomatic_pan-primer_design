#!/usr/bin/env python

import argparse
import subprocess
import pandas as pd
import numpy as np
from os import path
from collections import Counter
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
from datetime import datetime
import shutil
import math
import matplotlib.pyplot as plt

DT_STRING = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def parse_args():
    parser = argparse.ArgumentParser(description='Create MAFFT alignment with entropy plot')
    parser.add_argument('-f','--fasta', help='Fasta file for alignment', required=True)
    parser.add_argument('-m','--method', help='Mafft method.(default: %(default)s)', default='auto', choices=['auto', 'einsi', 'linsi', 'ginsi'])
    parser.add_argument('-p', '--prefix', help='Prefix. Default is input fasta')
    parser.add_argument('-d', '--out-dir', help='Output directory. Default is current working directory')
    parser.add_argument('--aa', help='Amino acid mode', default=False, action='store_true')
    parser.add_argument('-z', '--threads', help='Number of CPU threads to use. (default: %(default)s)', default='20', type=str)
    parser.add_argument('-t', '--plot-title', help='Title of entropy plot. Default is input filename')
    parser.add_argument('--skip-aln', help='Skip alignmnet. -f is assumed to be aligned file.', default=False, action='store_true')
    parser.add_argument('-i', '--plot-interactive', help='Interactive plot (default: %(default)s)', default=False, action='store_true')
    parser.add_argument('--memsave', help='Use memsave instead of nomemsave when running MAFFT (default: %(default)s)', default=False, action='store_true') 

    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    args = parser.parse_args()
    return args


def runProcess(cmd, pipe=None):
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

    if pipe is None:
        process = subprocess.run(cmd, stderr = suppressLog())
    else:
        file_pipe = open(pipe, 'w')
        process = subprocess.run(cmd, stdout=file_pipe)
    if process.returncode != 0:
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
        prefix = path.basename(args.fasta)
        prefix = prefix.replace(f'.{prefix.split(".")[-1]}', '')
        if args.skip_aln:
            args.prefix = prefix
            args.method = None
        else:
            args.prefix = f'{prefix}_mafft-{args.method}'
            if args.plot_title is None:
                args.plot_title = prefix
    if args.plot_title is None:
        args.plot_title = args.prefix
    if args.out_dir is None:
        file_dir= os.getcwd()
    else:
        file_dir = path.abspath(args.out_dir)
    args.out_dir = f'{file_dir}/mafft-entropy_{args.prefix}_{DT_STRING}'
    args.post_dir = f'{args.out_dir}/post_processing'

    if not path.exists(args.fasta):
        sys.exit(f'{args.fasta} not found')
    args.fasta = path.abspath(args.fasta)

    args.threads = validateThreads(args.threads)

    scripts_path = os.path.dirname(os.path.realpath(__file__))
    args.rmd = f'{scripts_path}/mafft_entropy.Rmd'
    args.rscript = f'{scripts_path}/mafft_entropy.R'
    global RUN_LOG
    RUN_LOG = f'{scripts_path}/logs/mafft-entropy_{DT_STRING}.log'
    if not path.exists(path.dirname(RUN_LOG)):
        os.mkdir(path.dirname(RUN_LOG))

    return args

def createOutputDir(args):
    dir_list = [args.out_dir, args.post_dir]
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
    tab = pd.DataFrame(args_tuple)
    tab.columns = ['Argument', 'Value']
    output_config = args.prefix+'_mafft-entropy.config'
    
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


################# MAFFT #################
def mafftSTD(args):
    cmd = ['mafft','--auto', '--reorder', '--anysymbol', '--adjustdirection', '--thread', args.threads]
    if args.memsave:
        cmd.extend(['--memsave'])
    else:
        cmd.extend(['--nomemsave'])
    if args.aa:
        cmd.extend(['--amino'])
    else:
        cmd.extend(['--nuc'])
    cmd.extend([args.fasta])
    aln_file = f'{args.out_dir}/{args.prefix}.aln'
    runProcess(cmd, pipe=aln_file)
    return aln_file

def mafftEGLinsi(args):
    num_records = sum(1 for _ in SeqIO.parse(args.fasta, 'fasta'))
    if args.method == 'einsi' and num_records>200:
        print(f'[WARNING] Computational time and memory usage can become significant, and alignment accuracy may decrease if the {num_records} are highly divergent. Up to ~200 highly divergent sequences can typically be handled effectively by E-INS-i')
    cmd = [args.method,'--reorder', '--anysymbol', '--adjustdirection', '--thread', args.threads]
    if args.memsave:
        cmd.extend(['--memsave'])
    else:
        cmd.extend(['--nomemsave'])
    if args.aa:
        cmd.extend(['--amino'])
    else:
        cmd.extend(['--nuc'])
    cmd.extend([args.fasta])
    aln_file = f'{args.prefix}.aln'
    runProcess(cmd, pipe=aln_file)
    return aln_file

def runMafft(args):
    if args.method == 'auto':
        aln_file = mafftSTD(args)
    else:
        aln_file =mafftEGLinsi(args)
    printInfo(f'Created {aln_file}')
    return aln_file

################# Shannon entropy #################
def trimEnds(args, aln_file):
    aln = AlignIO.read(aln_file, 'fasta')
    start = None
    end = None
    aln_len = aln.get_alignment_length()
    for i in range(0, aln_len):
        if all(rec.seq[i] !='-' for rec in aln):
            start = i
            break
    for i in range(aln_len-1, -1, -1):
        if all(rec.seq[i] !='-' for rec in aln):
            end = i+1
            break
    
    aln_trim = aln[:, start:end]
    output_file = f'{args.post_dir}/{args.prefix}_trim.aln'
    AlignIO.write(aln_trim, output_file, 'fasta')
    return aln_trim

def preprocessAln(args, aln_file):
    aln_tuple = []
    aln_trim = trimEnds(args, aln_file)
    for rec in aln_trim:
        seq = rec.seq.lower()
        if 'u' in seq:
            seq = seq.back_transcribe()
        aln_tuple.append((rec.id, str(seq)))
    return aln_tuple 

def calcEntropy(chars_list, states):
    ent = 0
    chars_len = len(chars_list)
    if chars_len <= 1:
        return ent
    _, counts = np.unique(chars_list, return_counts=True) #calculate the number of unique characters and their counts
    probabilties = counts/chars_len
    if np.count_nonzero(probabilties) <= 1:
        return ent
    for prob in probabilties:
        ent -= prob*math.log(prob, states)
    return ent

def alignmentEntropy(args,aln):
    printInfo('Calculating entropy')
    ent_pos_tuple = []
    for nuc_pos in range(0, len(aln[0][1])):
        pos = []
        for seq_num in range(0, len(aln)):
            name, seq = aln[seq_num]
            pos.append(seq[nuc_pos])
        if args.aa:
            ent = calcEntropy(pos, 20) #20 stats because sequence is aa
        else:
            ent = calcEntropy(pos, 4) #4 stats because sequence is nt
        ent_pos_tuple.append((nuc_pos, ent))
    ent_df = pd.DataFrame(ent_pos_tuple)
    ent_df.columns = ['pos0', 'ent']
    ent_df['mean'] = ent_df['ent'].rolling(10, center=True).mean()
    df_out = f'{args.post_dir}/{args.prefix}_trim_entropy.tsv'
    ent_df.to_csv(df_out, sep='\t', index=False)
    print(ent_df)
    return ent_df

def entropyPlot(args, aln):
    """
    Creates a single entropy plot.
    """
    fig_out = f'{args.out_dir}/{args.prefix}.svg'
    df1 = alignmentEntropy(args, aln)
    printInfo('Creating entropy plot')
    df1['pos1'] = df1['pos0']+1
    # Create the plot
    fig, ax = plt.subplots(figsize=(22, 4.8))  # Set figure size (WxH)

    ax.fill_between(df1["pos1"], df1["ent"], color="gainsboro", label="Entropy")
    ax.plot(df1["pos1"], df1["mean"], color="black", label="Rolling mean entropy", linewidth=0.5)
    ax.set_ylim((0, 1))
    ax.set_xlim(0, max(df1["pos1"]))
    ax.set_xlabel("Alignment position")
    ax.set_ylabel("Normalized Shannon's entropy")
    if args.method is None:
        ax.set_title(f"Mafft alignment\n{args.plot_title}")
    else:
        ax.set_title(f"Mafft {args.method} alignment\n{args.plot_title}")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.legend(loc=(0.835, 0.6))

    # Show the plot
    fig.savefig(fig_out, bbox_inches='tight')
    if args.plot_interactive:
        plt.show()
    else:
        plt.close()
    return fig_out


################# Main #################
def main():
    args = parse_args()
    args = checkArgs(args)
    printMessage('MAFFT alignment and entropy plot')
    createConfig(args, var_dict={}, print_config=True) 
    createOutputDir(args)
    if args.skip_aln:
        printInfo('Provided file is alignment file')
        aln_file = f'{args.out_dir}/{path.basename(args.fasta)}'
        shutil.copyfile(args.fasta, aln_file)
    else:
        aln_file = runMafft(args)
    aln_trim = preprocessAln(args, aln_file)
    plot = entropyPlot(args, aln_trim)
    _ = createConfig(args, var_dict={'mafft_aln':aln_file, 'plot':plot})
    print('Done')

if __name__ == "__main__":
   main()       