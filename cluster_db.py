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
import random
from sklearn.metrics import silhouette_score
import math
import polars as pl
from tqdm import tqdm

DT_STRING = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def parse_args():
    parser = argparse.ArgumentParser(description='Cluster database')
    parser.add_argument('-f','--fasta', help='Fasta file containing the original database', required=True)
    parser.add_argument('-m','--meta', help='TSV file containing meta-data for the fasta file, such as "species"', required=True)
    parser.add_argument('-c','--major-cluster', help='Range for major clustering. (default: %(default)s)', default='0.70 0.75 0.80 0.85', type=str)
    parser.add_argument('-s','--subcluster', help='Range for sub-clustering. (default: %(default)s)', default='0.96 0.97 0.98 0.99', type=str)
    parser.add_argument('-n','--numseq', help='Number of sequences per major cluster (default: %(default)s)', default=100, type=int)
    parser.add_argument('-p', '--prefix', help='Prefix. Default is input fasta')
    parser.add_argument('-o', '--out-dir', help='Output directory. Default is current working directory')
    parser.add_argument('-z', '--threads', help='Number of CPU threads to use. (default: %(default)s)', default='20')
    parser.add_argument('--max-len', help='Maximum sequence length for vsearch. (default: %(default)s)', default='50000')
    parser.add_argument('--all', help='Subcluster all major clusters regardless of silhouette score selection (default: %(default)s)', default=False, action='store_true')
    parser.add_argument('--overcluster-cutoff', help='Sets the cut-off for overclustering penalty for calculating the max acceptable sub_id (default: %(default)s)', default=0.2, type=int)
    parser.add_argument('--align-dir', help="Provide the directory of the vclust alignment to skip this step. This is useful if vclust align has already been performed on the input database")
    parser.add_argument('--mash', help="Provide a precomputed mash distance matrix file to skip the matrix generation step. This is useful if you've already run <mash dist -t> on the input database and want to reuse the output. (default: %(default)s)")
    parser.add_argument('--report-title', help='Title of report. Default is input filename')
    parser.add_argument('--skip-major', help='Skip major clustering by recreating the clusters from the provided cluster tsv')
    parser.add_argument('--cluster-method', help='Method to use for clustering (default: %(default)s)', default='vclust-single', type=str, choices=['vclust-single', 'vclust-complete', 'vclust-uclust', 'vclust-cd-hit', 'vclust-set-cover', 'vclust-leiden', 'vsearch', 'cd-hit-est'])
    parser.add_argument('-a', '--preferred-accessions', help="Path to file containing accessions to preferentially retain during database reduction. For each subcluster listed accessions are selected first (if present) and remaining slots are filled by random sampling ")

    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    # parser.add_argument('-', '--', help='')
    args = parser.parse_args()
    return args

def runProcess(cmd, pipe=None, suppress=False, print_cmd=True):
    if print_cmd:
        base_names = []
        for item in cmd:
            if item.endswith('/'):
                item_split = item.split('/')
                basename = f'{item_split[len(item_split)-2]}/'
            else:
                basename = path.basename(item)
            base_names.append(basename)
        # base_names= [path.basename(item) for item in cmd]
        print_cmd = f"\n\t{' '.join(base_names)}\n"
        print(print_cmd)
    log_out = suppressLog()
    log_out.write(f'\n{" ".join(cmd)}\n')
    log_out.close()

    if suppress:
        process = subprocess.run(cmd, stdout=suppressLog(), stderr = suppressLog())
    elif pipe is None:
        process = subprocess.run(cmd, stderr = suppressLog())
    else:
        file_pipe = open(pipe, 'w')
        process = subprocess.run(cmd, stdout=file_pipe)
    if process.returncode != 0:
        fail_msg = f'{cmd} FAILED'
        quitLog(fail_msg)


def roundHalfUp(n, decimals=0):
    if pd.isna(n):
        return np.nan
    multiplier = 10**decimals
    return math.floor(n * multiplier + 0.5) / multiplier

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

def printLog(msg):
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

def validateValues(var, name):
    values = [x.strip() for x in var.split()]
    try:
        floats = [float(x) for x in values]
    except:
        sys.exit(f'[ERROR] Non-numeric value found in {name}: {values}')
    
    out_of_range = [v for v in floats if not (0<=v<=1)]
    if out_of_range:
        sys.exit(f'[ERROR] {name} contains values outside of valid range [0, 1]: {out_of_range}')

    str_list = [str(v) for v in floats]
    return str_list

def findAlignFiles(args):
    required_endings = {'align.tsv': 'vclust_align', 'align.ids.tsv': 'vclust_ids'}
    files = os.listdir(args.align_dir)

    # Match and assign paths for each required ending
    for ending, attr in required_endings.items():
        matched = [f for f in files if f.endswith(ending)]
        if not matched:
            sys.exit(f"[ERROR] Cannot find file ending in '{ending}' in {args.align_dir}")
        elif len(matched) > 1:
            sys.exit(f"[ERROR] Multiple files end with '{ending}'")
        setattr(args, attr, os.path.join(args.align_dir, matched[0]))
    return args

def checkArgs(args):
    if args.prefix is None:
        prefix = path.basename(args.fasta)
        args.prefix = f'{path.splitext(prefix)[0]}_{args.cluster_method}'
    if args.report_title is None:
        args.report_title = args.prefix
    if args.out_dir is None:
        file_dir= os.getcwd()
    else:
        file_dir = path.abspath(args.out_dir)
    args.out_dir = f'{file_dir}/cluster-db_{args.prefix}_{DT_STRING}'
    if not path.exists(args.fasta):
        sys.exit(f'{args.fasta} not found')
    args.fasta = path.abspath(args.fasta)
    if not path.exists(args.meta):
        sys.exit(f'{args.meta} not found')
    if args.mash is not None:
        if not path.exists(args.mash):
            sys.exit(f'{args.mash} not found')
        args.mash = path.abspath(args.mash)
    if args.skip_major is not None:
        if not path.exists(args.skip_major):
            sys.exit(f'{args.skip_major} not found')
        args.skip_major = path.abspath(args.skip_major)
    if args.preferred_accessions is not None:
        if not path.exists(args.preferred_accessions):
            sys.exit(f'{args.preferred_accessions} not found')
        args.preferred_accessions = path.abspath(args.preferred_accessions)

    args.meta = path.abspath(args.meta)
    args.major_cluster = validateValues(args.major_cluster, '--major-cluster')
    args.subcluster = validateValues(args.subcluster, '--subcluster')

    if args.cluster_method == 'cd-hit-est':
        if any(float(x) < 0.80 for x in args.major_cluster):
            sys.exit(f'[ERROR] cd-hit-est does not function below 0.80')
    elif 'vclust' in args.cluster_method:
        if args.mash is not None:
            sys.exit(f'[ERROR] vclust clustering cannot be used with --mash')
    
    if args.align_dir is not None:
        args.align_dir = path.abspath(args.align_dir)    
        args = findAlignFiles(args)

    args.threads = validateThreads(args.threads)
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    args.rmd = f'{scripts_path}/cluster_db.Rmd'
    args.rscript = f'{scripts_path}/cluster_db.R'
    args.temp_dir = f'{args.out_dir}/temp'
    args.reduce_dir = f'{args.out_dir}/4_reduced_db'
    args.cluster_dir = f'{args.out_dir}/2_clustering'
    args.scores_dir = f'{args.out_dir}/3_scores'
    args.dist_dir = f'{args.out_dir}/1_dist'
    global RUN_LOG
    RUN_LOG = f'{scripts_path}/logs/cluster-db_{DT_STRING}.log'
    if not path.exists(path.dirname(RUN_LOG)):
        os.mkdir(path.dirname(RUN_LOG))
    return args

def createOutputDir(args):
    dir_list = [args.out_dir, args.reduce_dir, args.cluster_dir, args.scores_dir, args.dist_dir, args.temp_dir]
    for item in dir_list:
        os.makedirs(item)
    printInfo(f'\nOutput directory: {args.out_dir}')
    os.chdir(args.out_dir)

def readMeta(args):
    printInfo('Wrangling metadata')
    ### Validate that metadata file contains required columns."""
    required_options = [
        {'accession', 'species'},
        {'accession', 'viral_species'}
    ]

    # --- Load metadata ---
    meta_df = pd.read_csv(args.meta, sep='\t', header=0)
    # --- Check columns ---
    meta_cols = set(meta_df.columns)

    if not any(required.issubset(meta_cols) for required in required_options):
        quitLog(
            f"[ERROR] Metadata must contain either columns "
            f"{list(required_options[0])} or {list(required_options[1])}.\n"
            f"Found columns: {sorted(meta_cols)}"
        )
    meta_df = meta_df.dropna(subset=['accession'])
    cols_to_drop = {'cluster', 'major_id', 'sub_cluster', 'sub_id'}
    meta_df = meta_df.drop(columns=cols_to_drop & set(meta_df.columns))

    if 'viral_species' in meta_df.columns and 'species' not in meta_df.columns:
        meta_df = meta_df.rename(columns={'viral_species': 'species'})
    elif 'viral_species' in meta_df.columns and 'species' in meta_df.columns:
        meta_df = meta_df.rename(columns={'species':'viral_name'})
        meta_df = meta_df.rename(columns={'viral_species': 'species'})

    accessions_fasta = {rec.id for rec in SeqIO.parse(args.fasta, "fasta")}
    meta_df = meta_df[meta_df["accession"].isin(accessions_fasta)].copy()
    accessions_meta = set(meta_df['accession'].to_list())
    missing = sorted(accessions_fasta - accessions_meta)
    if accessions_fasta != accessions_meta:
        missing_df = pd.DataFrame({"missing_accession": missing})
        quitLog(f'[ERROR] All fasta records are not in meta_df. Missing:\n\n{missing_df}')
    return meta_df, accessions_fasta

def mashSketch(args, fasta, sketch_out):
    cmd = ['mash', 'sketch', '-i', fasta, '-o', sketch_out, '-k', '12', '-p', args.threads]
    runProcess(cmd)

def mashDist(args, mash_sketch, dist_out):
    cmd = ['mash', 'dist','-p', args.threads, '-t', mash_sketch, mash_sketch]
    runProcess(cmd, pipe=dist_out)

def runMash(args, fasta, sketch_out, dist_out):
    mashSketch(args, fasta, sketch_out)
    mashDist(args, sketch_out, dist_out)

def mashDB(args):
    #runs mash on the original database. If precomputed mash is provided, then that is used instead
    if args.mash is None:
        printMessage(f'Generating Mash dist matrix')
        sketch_out = f'{args.dist_dir}/{args.prefix}_sketch.msh'
        dist_out = f'{args.dist_dir}/{args.prefix}_mash.tsv'
        runMash(args, args.fasta, sketch_out, dist_out)
    else:
        dist_out = f'{args.dist_dir}/{args.prefix}_mash.tsv'
        printInfo(f'Copying mash dist file')
        shutil.copy2(args.mash, dist_out)
    
    printInfo('Reading mash dist matrix')
    mash_df = pl.read_csv(dist_out, separator="\t", has_header=True).rename({'#query':'query'})
    printLog(mash_df)
    return mash_df

def retrieveWordSize(id):
    id = float(id)
    if id >=0.975:
        word_size = 11
    elif id >=0.95:
        word_size = 10
    elif id >= 0.925:
        word_size = 9
    elif id >= 0.90:
        word_size = 8
    elif id >= 0.88:
        word_size = 7
    elif id >= 0.85:
        word_size = 6
    elif id >= 0.80:
        word_size = 5
    return str(word_size)

def cdHit(args, fasta, id, out_dir, rec_dict, print_cmd=False):
    #runs cd-hist-est on the provided fasta. Output is placed in temp dir, which is then removed after parsing
    out_prefix = f'{out_dir}/temp'
    word_size = retrieveWordSize(id)
    cmd = ['cd-hit-est', '-i', fasta, '-o', out_prefix, '-c', id, '-n', word_size, '-d', '0', '-T', args.threads, '-M', '0']
    runProcess(cmd, suppress=True, print_cmd=print_cmd)
    
    # print('   ...finding clusters...')
    cluster_dict = dict()
    with open(f'{out_prefix}.clstr', 'r') as fp:
        for line in fp:
            if line.startswith('>Cluster'):
                   cluster = int(line.replace('>Cluster', '').strip())
            else:
                acc = line.split('>')[1].split('...')[0]
                rec = rec_dict[acc]
                if cluster not in cluster_dict:
                    cluster_dict[cluster] =[rec]
                else:
                    cluster_dict[cluster].append(rec)
    cluster_tuple = []
    # print('   ...creating clustered records...')
    for key in cluster_dict.keys():
        records = cluster_dict[key]
        records_out = f'{out_dir}/{key}.fasta'
        SeqIO.write(records, records_out, 'fasta')
        cluster_tuple.append((key, records_out))
    os.remove(out_prefix)
    os.remove(f'{out_prefix}.clstr')
    return cluster_tuple

def retrieveRecords(args):
    rec_dict = dict()
    with open(args.fasta, "r") as gb_file:
        for rec in SeqIO.parse(gb_file, 'fasta'):
            acc = rec.id
            rec_dict[acc] = rec
    return rec_dict

def parseClusters(cluster_tuple, major_id, cluster_dir):
    major_tuple = []
    for cluster, fasta in cluster_tuple:
        with open(fasta) as fopen:
            for rec in SeqIO.parse(fopen, 'fasta'):
                acc = rec.id
                major_tuple.append((acc, major_id, cluster, fasta, cluster_dir))
            
    major_df = pd.DataFrame(major_tuple)
    major_df.columns=['accession', 'major_id', 'cluster', 'fasta', 'cluster_dir']
    return major_df

def runCluster(args, fasta, id, prefix_dir, rec_dict):
    if 'vclust' in args.cluster_method:
        cluster_dir = f'{prefix_dir}/{args.cluster_method}_{id}'
        os.makedirs(cluster_dir)
        cluster_df = vClustCluster(args, args.vclust_align, args.vclust_ids, id, cluster_dir)
    elif args.cluster_method == 'cd-hit-est':
        cluster_dir = f'{prefix_dir}/cdhit_{id}'
        os.makedirs(cluster_dir)
        cluster_tuple = cdHit(args, fasta, id, cluster_dir, rec_dict)
        cluster_df = parseClusters(cluster_tuple, id, cluster_dir)
    elif args.cluster_method == 'vsearch':
        cluster_dir = f'{prefix_dir}/vsearch_{id}'
        os.makedirs(cluster_dir)
        cluster_tuple = vSearch(args, fasta, id, cluster_dir)
        cluster_df = parseClusters(cluster_tuple, id, cluster_dir)

    return cluster_df

def vSearch(args, fasta, id, vsearch_dir, print_cmd=False):
    cmd = ['vsearch', '--cluster_fast', fasta, '--clusters', f'{vsearch_dir}/', '--id', id, '--threads', args.threads, '--maxseqlength', args.max_len]
    runProcess(cmd, print_cmd=print_cmd)
    cluster_files = [path.join(vsearch_dir, item) for item in os.listdir(vsearch_dir)]
    cluster_tuple = []
    for file in cluster_files:
        cluster = path.basename(file)
        new_name = f'{file}.fasta'
        os.rename(file, new_name)
        cluster_tuple.append((cluster, new_name))
    return cluster_tuple

def vClustAlign(args):
    printMessage('Vclust Average Nucleotide Identity (ANI)')
    vclust_align_file = f'{args.dist_dir}/{args.prefix}_align.tsv'
    vclust_ids_file = f'{args.dist_dir}/{args.prefix}_align.ids.tsv'
    if args.align_dir is None:
        printInfo(f'Generating all vs all ANI')
        printLog(f'   ...this may take a good while...')
        args.vclust_align = vclust_align_file
        args.vclust_ids = vclust_ids_file
        cmd = ['vclust', 'align', '-i', args.fasta, '-o', args.vclust_align]
        runProcess(cmd)
    else:
        printInfo('Copying files')
        shutil.copy2(args.vclust_align, vclust_align_file)
        shutil.copy2(args.vclust_ids, vclust_ids_file)
    printInfo('Reading ANI dataframe')
    ani_df = pl.read_csv(
            args.vclust_align,
            separator="\t",
            has_header=True,
            schema_overrides={
                "qidx": pl.UInt32,
                "ridx": pl.UInt32,
                "query": pl.String,
                "reference": pl.String,
                "tani": pl.Float32,
                "gani": pl.Float32,
                "ani": pl.Float32,
                "qcov": pl.Float32,
                "rcov": pl.Float32,
                "num_alns": pl.UInt16,
                "len_ratio": pl.Float32,
            },
    )
    ani_ids = pl.read_csv(args.vclust_ids, separator='\t', has_header=True)
    printInfo("Extrapolating distance matrix")
    pivot_df = ani_df.select(["query", "reference", "tani"]).with_columns((1 - pl.col("tani")).alias("dist"))
    
    mat = pivot_df.pivot(index="query", on="reference", values="dist").fill_null(0.0)
    # Extract accessions in row order
    row_order = mat["query"].to_list()
    # Get all numeric columns (distance columns)
    num_cols = [c for c in mat.columns if c != "query"]
    # Filter only the accessions that exist as columns (some might be missing)
    cols_in_order = [c for c in row_order if c in num_cols]
    # Reorder columns to match the row order
    mat = mat.select(["query"] + cols_in_order)
    printLog(mat)

    return args, ani_df, ani_ids, mat
    
def vClustCluster(args, align_tsv, ids_tsv, id, vclust_dir, print_cmd=False):
    method = args.cluster_method.replace('vclust-','')
    cluster_file =  f'{vclust_dir}/vclust_{id}.tsv'
    cmd = ['vclust', 'cluster', '-i', align_tsv, '-o' , cluster_file, '--ids', ids_tsv , '--algorithm', method, '--metric', 'tani', '--tani', id]
    runProcess(cmd, print_cmd=print_cmd)
    df1 = pd.read_csv(cluster_file, sep='\t')
    df1.columns = ['accession', 'cluster']
    return df1

def subclusterAni(cluster_df, sub_id, ani_df, ani_ids, sub_dir, cluster):
    ani_out = f'{sub_dir}/{sub_id}_{cluster}_ani.tsv'
    ids_out = f'{sub_dir}/{sub_id}_{cluster}_ani.ids.tsv'
    df1 = cluster_df[cluster_df['cluster']==cluster]
    accessions = set(df1['accession'].to_list())
    ani_sub = ani_df.filter(pl.col("query").is_in(accessions) & pl.col("reference").is_in(accessions)).drop(["qidx", "ridx"], strict=False)
    ids_sub = ani_ids.filter(pl.col("id").is_in(accessions)).with_row_index("index")
    ani_sub = ani_sub.join(ids_sub.select(["id", "index"]).rename({"id": "query", "index": "qidx"}), on="query", how="left")
    ani_sub = ani_sub.join(ids_sub.select(["id", "index"]).rename({"id": "reference", "index": "ridx"}), on="reference", how="left")
    ani_sub = ani_sub.select(['qidx', 'ridx', 'query', 'reference', 'tani', 'gani', 'ani', 'qcov', 'rcov', 'num_alns', 'len_ratio'])
    ids_sub = ids_sub.select(["id", "seq_len", "no_parts"])
    ani_sub.write_csv(ani_out, separator="\t")
    ids_sub.write_csv(ids_out, separator="\t")

    return ani_out, ids_out

def subclusterVclust(args, cluster_df, major_id, ani_df, ani_ids, meta_df):
    subcluster_list = []
    sub_ids = args.subcluster
    clusters = sorted(set(cluster_df['cluster'].to_list()))
    total_steps = len(sub_ids)*len(clusters)
    desc=f'Subclustering major id {major_id}'

    with tqdm(total=total_steps, desc=desc, position=1, leave=True) as pbar:
        for sub_id in sub_ids:
            sub_dir =  f'{args.out_dir}/temp/{args.cluster_method}_{major_id}'
            if not path.exists(sub_dir):
                os.makedirs(sub_dir)
            for cluster in clusters:
                ani_out, ids_out = subclusterAni(cluster_df, sub_id, ani_df, ani_ids, sub_dir, cluster)
                df1 = vClustCluster(args, ani_out, ids_out, sub_id, sub_dir)
                df1 = df1.rename(columns={'cluster':'sub_cluster'})
                df1['cluster']= cluster
                df1["sub_cluster"] = str(cluster) + "-" + df1["sub_cluster"].astype(str)
                df1['sub_id'] = sub_id
                df1['major_id']=major_id
                subcluster_list.append(df1)

                pbar.update(1)  

    subcluster_df = pd.concat(subcluster_list, ignore_index=True)
    subcluster_df = subcluster_df[['accession', 'major_id', 'sub_id', 'cluster', 'sub_cluster']]
    subcluster_df= subcluster_df.merge(meta_df, on='accession', how='left')
    return subcluster_df

def subClusterCdhitVsearch(args, cluster_dir, cluster_dict, major_id, rec_dict, meta_df):
    sub_ids = args.subcluster
    clusters = sorted(list(cluster_dict.keys()))
    total_steps = len(sub_ids)*len(clusters)
    desc=f'Subclustering major id {major_id}'

    sub_df_list = []
    with tqdm(total=total_steps, desc=desc, position=1, leave=True) as pbar:
        for sub_id in args.subcluster:
            i = 0
            sub_dir =  f'{cluster_dir}/{args.cluster_method}_{sub_id}'
            os.makedirs(sub_dir)
            for cluster in clusters:
                fasta = cluster_dict[cluster]
                i += 1
                temp_dir = f'{sub_dir}/temp_{i}'
                os.makedirs(temp_dir)
                # cluster_tuple = cdHit(args, fasta, sub_id, temp_dir)
                if args.cluster_method == 'cd-hit-est':
                    _ = cdHit(args, fasta, sub_id, temp_dir,rec_dict)
                else:
                    _ = vSearch(args, fasta, sub_id, temp_dir)
                sub_files = [path.join(temp_dir, item) for item in os.listdir(temp_dir)]   
                sub_tuple = []
                for file in sub_files:
                    if file.endswith('.fasta'):
                        sub_cluster = f'{cluster}-{path.basename(file).split(".")[0]}'
                        with open(file) as fopen:
                            for rec in SeqIO.parse(fopen, 'fasta'):
                                acc = rec.id 
                                sub_tuple.append((acc, major_id, sub_id, cluster, sub_cluster))
                
                sub_df1 = pd.DataFrame(sub_tuple)
                sub_df1.columns = ['accession', 'major_id', 'sub_id', 'cluster', 'sub_cluster']
                sub_df_list.append(sub_df1)
                shutil.rmtree(temp_dir)
                pbar.update(1)  
            shutil.rmtree(sub_dir)
    sub_df2 = pd.concat(sub_df_list, ignore_index=True)
    sub_df2 = sub_df2.merge(meta_df, on='accession', how='left')
    return sub_df2

def skipClustering(args, rec_dict):
    printInfo('Reading clusters')
    df1 = pd.read_csv(args.skip_major, sep='\t', header=0)
    major_ids = set(df1['major_id'].to_list())
    args_ids = set(map(float, args.major_cluster))
    missing_ids = [x for x in args_ids if x not in major_ids]
    common_ids = sorted(major_ids & args_ids)
    if missing_ids:
        quitLog(f'[ERROR] {sorted(missing_ids)} not found in major cluster in provided tsv {sorted(major_ids)}')
    elif major_ids != args_ids:
        printWarning(f'Not the same major clusters in {sorted(args_ids)} and provided tsv {sorted(major_ids)}. Will continue with {common_ids}')
    df1['major_id'] = df1['major_id'].astype(str)
    df1['cluster'] = df1['cluster'].astype(str) 

    if 'vclust' in args.cluster_method:
        df2 = df1[df1['accession', 'major_id', 'cluster']].copy()
    else:
        printLog(f'   ...recreating clustered records...')
        fasta_tuple = []
        for major_id in sorted(set(map(str, common_ids))):
            printLog(f'   ...{major_id}...')
            cluster_dir = f'{args.temp_dir}/cluster_{major_id}'
            os.makedirs(cluster_dir)
            major_df = df1[df1['major_id']==major_id].copy()
            clusters = sorted(set(major_df['cluster'].to_list()))
            for cluster in clusters:
                cluster_df = major_df[major_df['cluster']==cluster].copy()
                acc_set = sorted(set(cluster_df['accession'].to_list()))
                cluster_rec = []
                fasta = f'{cluster_dir}/{cluster}.fasta'
                for acc in acc_set:
                    cluster_rec.append(rec_dict[acc])
                    fasta_tuple.append((acc, major_id, cluster, fasta, cluster_dir))
                SeqIO.write(cluster_rec, fasta, 'fasta')
        
        fasta_df = pd.DataFrame(fasta_tuple)
        fasta_df.columns = ['accession', 'major_id', 'cluster', 'fasta', 'cluster_dir']
        df2 = df1.merge(fasta_df, on=['accession', 'major_id', 'cluster'], how='left')
        df2 = df2[['accession', 'major_id', 'cluster', 'fasta', 'cluster_dir']]
    return df2
    
def majorClusters(args, meta_df):
    printMessage('Clustering with major cluster ids')
    rec_dict = retrieveRecords(args)
    if args.skip_major is None:
        major_tuple = []
        major_ids = args.major_cluster
        for major_id in tqdm(major_ids, desc=f'Major clusters {args.cluster_method}', position=0, leave=True):
            # printInfo(f'Clustering with {major_id}')
            df1 = runCluster(args, args.fasta, major_id, args.temp_dir, rec_dict)
            df1['major_id']=major_id
            major_tuple.append(df1)
        printMessage('Collating major clusters')                
        major_df = pd.concat(major_tuple, ignore_index=True)
        # major_df = major_df[['accession', 'major_id', 'cluster']]
    else:
        major_df = skipClustering(args, rec_dict)

    major_df = major_df.merge(meta_df, on='accession', how='left')
    major_out_file = f'{args.cluster_dir}/{args.prefix}_major-cluster.tsv'
    if 'fasta' in major_df.columns:
        major_out = major_df.copy().drop(columns=['fasta', 'cluster_dir'])
    else:
        major_out = major_df.copy()
    major_out.to_csv(major_out_file, sep = '\t', index=False)
    printLog(major_out)
    return major_df, rec_dict, major_out_file

def calculateSilhouette(dist_df, df1, col_name):
    dist_acc = set(dist_df['query'].to_list()) #which accession in dist
    cluster_acc = set(df1['accession'].to_list()) #which accessions in cluster
    not_found_dist = sorted(dist_acc-cluster_acc)
    not_found_cluster = sorted(cluster_acc-dist_acc)
    if len(not_found_dist) !=0:
        quitLog(f'[ERROR] Inconsistencies between DISTANCE MATRIX and clustering!\nThe following sequences are found in matrix but not in the clustering: {not_found_dist}\n')
    if len(not_found_cluster) !=0:
        quitLog(f'[ERROR] Inconsistencies between DISTANCE MATRIX and clustering!\nThe following sequences are found in clustering but not in the matrix: {not_found_cluster}\n')
    cluster_dict = dict(zip(df1['accession'], df1[col_name]))
    dist_df = dist_df.with_columns(pl.col("query").replace(cluster_dict).alias("cluster"))
    clusters = dist_df['cluster'].to_numpy()
    dist_matrix = dist_df.drop(["query", "cluster"]).to_numpy()
    clust_score = silhouette_score(dist_matrix, clusters, metric='precomputed')
    return clust_score

def majorScore(args, dist_df, major_df, meta_df):
    printMessage("Calculating silhouette score for major clusters")
    score_tuple = []
    for major_id in sorted(set(major_df['major_id'].to_list())):
        df1 = major_df[major_df['major_id']==major_id].copy()
        n_clusters = len(set(df1['cluster'].to_list()))
        if n_clusters < 2:
            printWarning(f'Cannot calculate score for {major_id} as there are only {n_clusters} clusters')
            score = np.nan
        else:
            score = calculateSilhouette(dist_df, df1, col_name='cluster')
        
        score_tuple.append((major_id, score))
    
    # calculate silhouette score based on species demarcation
    species_score = calculateSilhouette(dist_df, meta_df, col_name='species')
    score_tuple.append(('species', species_score))

    scores_df = pd.DataFrame(score_tuple)
    scores_df.columns = ['major_id', 'score']
    scores_df['score'] = scores_df['score'].apply(lambda x: roundHalfUp(x, 3))
    scores_out = f'{args.scores_dir}/{args.prefix}_major-silhouette-scores.tsv'
    scores_df.to_csv(scores_out, sep='\t', index=False)
    printLog(scores_df)
    scores_df = scores_df[scores_df['major_id'] != 'species']
    return scores_df, scores_out

def bestScoreSub(args, score_df, best_major):
    printInfo(f'Finding best subclustering for {best_major}')
    df1 = score_df[score_df['major_id']==float(best_major)].copy()
    df1['sub_id'] = df1['sub_id'].astype(float)
    
    if (df1['w_dev'] == 0).any():
        df2 = df1[df1['w_dev'] == 0].copy()
        max_score = df2['sub_id'].max() # want the highest subscluster id as w_dev is 0
    else:
        min_val = df1['w_dev'].min()
        df2 = df1[df1['w_dev'] == min_val].copy()
        max_score = df2['sub_id'].min() # want the lowest subcluster id as w_dev is above 0
       
    printLog(f'   ...best is {max_score}...')
    accept_df = df1[df1['w_dev']<=args.overcluster_cutoff].copy()
    accept_score = accept_df['sub_id'].max()
    if not (df1['w_dev'] == 0).any():
        printLog(f'   [WARNING] {max_score} is higher than 0. Might want to run with a range of lower subcluster ids')
    elif accept_score != max_score:
        printLog(f'   ...however, {accept_score} might be acceptable...')

    return max_score, accept_score

def bestMajor(major_scores_df, major_df):
    printInfo('Finding highest major score')
    max_score = major_scores_df['score'].max()
    best_majors = major_scores_df.loc[major_scores_df['score'] == max_score, 'major_id'].to_list()
    if len(best_majors) == 1:
        best_major = best_majors[0]
    else:
        printWarning(f'Multiple high scores found: {best_majors}')
        floats = []
        for ele in best_majors:
            floats.append(float(ele))
        best_major = str(max(floats))
    printLog(f'   ...score is {max_score} at major cluster {best_major}...')
    best_df = major_df[major_df['major_id']==best_major].copy()

    return best_major, best_df

def bestMajorSubVclust(args, best_df, best_major, ani_df, ani_ids, meta_df):
    sub_df = subclusterVclust(args, best_df, best_major, ani_df, ani_ids, meta_df)
    printLog(sub_df)
    sub_out = f'{args.cluster_dir}/{args.prefix}_sub-cluster.tsv'
    sub_df.to_csv(sub_out, sep='\t', index=False)

    return sub_df, sub_out

def bestMajorSubCdhitVsearch(args, best_df, best_major, rec_dict, meta_df):
    best_dict = dict(zip(best_df['cluster'], best_df['fasta']))
    cluster_dir = best_df['cluster_dir'].to_list()[0]
    printMessage('Subclustering')
    sub_df = subClusterCdhitVsearch(args, cluster_dir, best_dict, best_major, rec_dict, meta_df)
    printLog(sub_df)
    sub_out = f'{args.cluster_dir}/{args.prefix}_sub-cluster.tsv'
    sub_df.to_csv(sub_out, sep='\t', index=False)

    return sub_df, sub_out

def subSilhouetteScores(args, best_df, best_major, sub_df):
    printInfo('Calculating score for subclusters')
    scores_list = []
    for cluster in sorted(set(best_df['cluster'].to_list())):
        score_df = scoreSub(args, best_major, cluster, sub_df)
        scores_list.append(score_df)
    
    scores_df = pd.concat(scores_list, ignore_index=True)

    printInfo('Calculating weighted scores')
    weighted_scores = weightedSubScore(scores_df)
    weighted_scores['major_id'] = weighted_scores['major_id'].astype(float)
    weighted_scores['sub_id'] = weighted_scores['sub_id'].astype(float)
    weighted_scores = weighted_scores.sort_values(by=['major_id', 'sub_id'], ignore_index=True)
    scores_out = f'{args.scores_dir}/{args.prefix}_sub-scores.tsv'
    weighted_scores.to_csv(scores_out, sep='\t', index=False)
    printLog(weighted_scores)
    return weighted_scores, scores_out

def allSubVclust(args, major_df, ani_df, ani_ids, meta_df):
    printMessage(f'Subclustering on all major ids')
    major_ids = sorted(set(major_df['major_id'].to_list()))
    sub_list = []
    for major_id in tqdm(major_ids, desc='Major clusters', position=0, leave=True):
        df1 = major_df[major_df['major_id']==major_id].copy()
        df2 = subclusterVclust(args, df1, major_id, ani_df, ani_ids, meta_df)
        sub_list.append(df2)

    sub_df = pd.concat(sub_list, ignore_index=True)
    sub_out = f'{args.cluster_dir}/{args.prefix}_sub-cluster.tsv'
    sub_df.to_csv(sub_out, sep='\t', index=False)

    return sub_df, sub_out


def allSubScores(args, sub_df):
    printInfo('Calculating weighted scores')
    weighted_list = []
    major_ids = sorted(set(sub_df['major_id'].to_list()))
    for major_id in major_ids:
        df1 = sub_df[sub_df['major_id']==major_id]
        clusters = sorted(set(df1['cluster'].to_list()))
        scores_list = []
        for cluster in clusters:
            score_df = scoreSub(args, major_id, cluster, df1)
            scores_list.append(score_df)

        scores_df = pd.concat(scores_list, ignore_index = True)
        weighted_scores = weightedSubScore(scores_df)
        weighted_list.append(weighted_scores)
    
    weighted_df = pd.concat(weighted_list)
    weighted_df['major_id'] = weighted_df['major_id'].astype(float)
    weighted_df['sub_id'] = weighted_df['sub_id'].astype(float)
    weighted_df = weighted_df.sort_values(by=['major_id', 'sub_id'], ignore_index=True)
    scores_out = f'{args.scores_dir}/{args.prefix}_sub-scores.tsv'
    weighted_df.to_csv(scores_out, sep='\t', index=False)

    printLog(weighted_df)
    return weighted_df, scores_out   


def allSubCdhitVsearch(args, major_df, rec_dict, meta_df):
    sub_list = []
    major_ids = sorted(set(major_df['major_id'].to_list()))
    for major_id in tqdm(major_ids, desc='Major clusters', position=0, leave=True):
        df1 = major_df[major_df['major_id']==major_id].copy()
        cluster_dict = dict(zip(df1['cluster'], df1['fasta']))
        cluster_dir = df1['cluster_dir'].to_list()[0]
        sub_df = subClusterCdhitVsearch(args, cluster_dir, cluster_dict, major_id, rec_dict, meta_df)
        sub_list.append(sub_df)

    sub_dfs = pd.concat(sub_list)
    sub_out = f'{args.cluster_dir}/{args.prefix}_sub-cluster.tsv'
    sub_dfs.to_csv(sub_out, sep='\t', index=False)

    return sub_dfs,sub_out
    
def scoreSub(args, major_id, cluster, sub_df):
    major_df = sub_df[sub_df['major_id']==major_id]
    cluster_df = major_df[major_df['cluster']==cluster]
    score_tuple = []
    for sub_id in sorted(set(cluster_df['sub_id'].to_list())):
        # print(f'...{sub_id}')
        df1 = cluster_df[cluster_df['sub_id']==sub_id]
        n_sequences = len(df1['cluster'].to_list())
        n_subclusters = len(set(df1['sub_cluster'].to_list()))
        if n_subclusters < args.numseq:
            deviation = 0
        else:
            deviation = abs(n_subclusters-args.numseq)/args.numseq

        score_tuple.append((major_id, sub_id, cluster, deviation, n_sequences, n_subclusters))
   

    scores_df = pd.DataFrame(score_tuple)
    scores_df.columns = ['major_id','sub_id', 'cluster', 'dev', 'n_sequences', 'n_subclusters']
    return scores_df

def weightedSubScore(df1):
    sub_ids = set(df1['sub_id'].to_list())
    major_ids = set(df1['major_id'].to_list())
    weighted_scores = []
    for major_id in major_ids:
        major_df = df1[df1['major_id']==major_id]
        for sub_id in sub_ids:
            sub_df = major_df[major_df['sub_id']==sub_id]
            total_weighted_dev = (sub_df['dev'] * sub_df['n_sequences']).sum() / sub_df['n_sequences'].sum()
            mean_dev = sub_df['dev'].mean()
            weighted_scores.append((major_id, sub_id, total_weighted_dev, mean_dev))
    
    weighted_df = pd.DataFrame(weighted_scores)
    weighted_df.columns = ['major_id', 'sub_id', 'w_dev', 'mean_dev']
    weighted_df['w_dev'] = weighted_df['w_dev'].apply(lambda x: roundHalfUp(x, 3))
    weighted_df['mean_dev'] = weighted_df['mean_dev'].apply(lambda x: roundHalfUp(x, 3))
    return weighted_df

def pickSeqs(args, subcluster_tuple):
    max_pick = args.numseq
    pick_dict = dict()  # Final selected clusters and sequences
    remains = []  # Remaining clusters and sequences after each iteration
    # Step 1: Ensure each cluster contributes at least one sequence
    for cluster, val in subcluster_tuple:
        pick_dict[cluster] = 1 # Each cluster contributes 1 sequence
        # Store the remaining sequences for further selection
        if val > 1:
            remains.append((cluster, val - 1))
    # Remaining sequences to pick after ensuring 1 per cluster
    remaining_total = max_pick - len(pick_dict)
    # Step 2-4: Continue picking until reaching the required number of sequences
    while remaining_total > 0 and remains:
        next_remains = []
        # Step 2: Go through all remaining clusters once
        avg_sequences = remaining_total / len(remains)  # Recalculate after each full loop
        for cluster, val in remains:
            if val <= avg_sequences:
                pick_dict[cluster] += val
                remaining_total -= val
            else:
                # Temporarily pick the average
                pick_amount = int(avg_sequences)
                if pick_amount == 0:
                    pick_amount = 1
                pick_dict[cluster] += pick_amount
                remaining_total -= pick_amount
                # Save the remaining sequences for the next round
                if val - pick_amount > 0:
                    next_remains.append((cluster, val - pick_amount))
        # Step 3: Update remains after a full loop
        remains = next_remains
    # df = pd.DataFrame(list(pick_dict.items()), columns=['Cluster', 'Picks'])
    # print(df)
    return pick_dict

def randomSeq(args, subcluster_df, preferred=None):
    cluster_sum = subcluster_df['cluster'].value_counts().reset_index()
    random_list = []
    for cluster in sorted(set(cluster_sum['cluster'].to_list())):
        df1 = subcluster_df[subcluster_df['cluster']==cluster][['accession','sub_cluster']]
        df1_sum = df1['sub_cluster'].value_counts().reset_index()
        df1_tuple = list(df1_sum.itertuples(index=False, name=None))
        pick_dict = pickSeqs(args, df1_tuple)
        for subcluster in pick_dict.keys():
            acc_list = subcluster_df[subcluster_df['sub_cluster']==subcluster]['accession'].to_list()
            num = pick_dict[subcluster]
            if preferred:
                acc_set = set(acc_list)
                forced_accessions = preferred.intersection(acc_set)
                random_list.extend(list(forced_accessions))
                num = num-len(forced_accessions)
                if num <1:
                    continue
                acc_list = list(acc_set-forced_accessions)
            random_accs = random.sample(acc_list, num)
            random_list.extend(random_accs)
        random_df = subcluster_df[subcluster_df['accession'].isin(random_list)]
    return random_df

def loadPreferred(args, accessions_fasta):
    """
    Read accession IDs from file (one accession per line).
    Blank lines are ignored.
    Raises ValueError if a line contains multiple fields or whitespace.
    """
    accessions = set()
    with open(args.preferred_accessions, "r") as f:
        for lineno, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            # Ignore blank lines
            if not line:
                continue
            # Reject whitespace-separated multiple entries
            if len(line.split()) != 1:
                quitLog(f"[ERROR] Invalid format in {path} (line {lineno}). File must contain exactly one accession per line.")
            if line not in accessions_fasta:
                quitLog(f'[ERROR] {line} from preferred accessions not present in fasta file')
            accessions.add(line)
 
    return accessions

def reduceDB(args, meta_df, df1, preferred=None):
    for major_id in set(df1['major_id'].to_list()):
        major_df = df1[df1['major_id']==major_id]
        for sub_id in set(major_df['sub_id'].to_list()):
            sub_df = major_df[major_df['sub_id']==sub_id]
            random_df = randomSeq(args, sub_df, preferred=preferred)
            # random_meta = random_df.merge(meta_df, on='accession', how='left')
            random_out = f'{args.reduce_dir}/{args.prefix}_{major_id}_{sub_id}_reduced.tsv'
            # random_meta.to_csv(random_out, sep='\t', index=False)
            random_df.to_csv(random_out, sep='\t', index=False)

    
def rmClusterFiles(args):
    printInfo(f'Removing {args.temp_dir}')
    shutil.rmtree(args.temp_dir)


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
    output_config = args.prefix+'_cluster-db.config'
    
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

def main():
    args = parse_args()
    args = checkArgs(args)
    createConfig(args, var_dict={}, print_config=True) 
    createOutputDir(args)
    meta_df, accessions_fasta = readMeta(args)
    if args.preferred_accessions:
        preferred_accessions = loadPreferred(args, accessions_fasta)
    if 'vclust' in args.cluster_method:
        args, ani_df, ani_ids, dist_df = vClustAlign(args)
    else:
        dist_df = mashDB(args)

    major_df, rec_dict, major_tsv= majorClusters(args, meta_df)
    major_scores_df, major_scores_tsv = majorScore(args, dist_df, major_df, meta_df)
    var_dict = add2varDict(['major_tsv', 'major_scores'], [major_tsv, major_scores_tsv])
    if args.all:
        if 'vclust' in args.cluster_method:
            sub_df, sub_out = allSubVclust(args, major_df, ani_df, ani_ids, meta_df)
        else:
            sub_df, sub_out = allSubCdhitVsearch(args, major_df, rec_dict, meta_df)
        weighted_scores, scores_out =  allSubScores(args, sub_df)
        var_dict = add2varDict(['sub_tsv', 'sub_scores'], [sub_out, scores_out], var_dict)
        
    else:
        best_major, best_df = bestMajor(major_scores_df, major_df)
        if 'vclust' in args.cluster_method:
            sub_df, sub_out = bestMajorSubVclust(args, best_df, best_major, ani_df, ani_ids, meta_df)
        else:
            sub_df, sub_out = bestMajorSubCdhitVsearch(args, best_df, best_major, rec_dict, meta_df)
        weighted_scores, scores_out = subSilhouetteScores(args, best_df, best_major, sub_df)
        max_sub, accept_sub = bestScoreSub(args, weighted_scores, best_major)
        var_dict = add2varDict(['sub_tsv', 'sub_scores', 'best_major', 'max_sub', 'accept_sub'], [sub_out, scores_out, best_major, max_sub, accept_sub], var_dict)

    if args.preferred_accessions:
        reduceDB(args, meta_df, sub_df, preferred=preferred_accessions)
    else:
        reduceDB(args, meta_df, sub_df)

    rmClusterFiles(args)
    output_config = createConfig(args, var_dict=var_dict)
    executeR(args, output_config)

    printLog('Done')
if __name__ == "__main__":
   main()       