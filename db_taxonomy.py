#!/usr/bin/env python

import pandas as pd
import argparse
from os import path
import sys
from Bio import SeqIO

KT_TAXONOMY_DB = '<USER-PATH>'

def parse_args():
    parser = argparse.ArgumentParser(description='Create taxonomy for blast db')
    parser.add_argument('-f','--fasta', help='Fasta file containing same accession as in meta', required=True)
    parser.add_argument('-m','--meta', help='Meta file containing acession and species', required=True)
    parser.add_argument('-c','--column', help='Meta data column for species. Default species', default='species')
    parser.add_argument('-o','--output', help='Output name (default <meta-file>_taxid.tsv')
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

def printInfo(message):
    msg = f'\n[INFO] {message}'
    print(msg)

def checkArgs(args):
    if args.output is None:
        output = path.basename(args.fasta)
        output = output.replace(f'.{output.split('.')[-1]}', '')
        args.output = f'{output}_taxid.tsv'
    if not path.exists(args.meta):
        sys.exit(f'Primer tsv file {args.meta} not found')
    if not path.exists(args.fasta):
        sys.exit(f'Primer tsv file {args.fasta} not found')
    args.meta = path.abspath(args.meta)
    args.fasta = path.abspath(args.fasta)
    
    return args

def taxonDict():
    printInfo('Preparing taxonomies')
    taxon_db = pd.read_csv(KT_TAXONOMY_DB, sep = '\t', header=0, names=['taxid', 'node', 'parent', 'level', 'name'])
    taxon_db = taxon_db[taxon_db['level']=='species']
    taxon_dict = taxon_db.set_index('name')['taxid'].to_dict()
    taxon_dict['unclassified']=0
    print('  Finished preparations')
    return taxon_dict


def dbTaxonomy(args):
    meta_data = pd.read_csv(args.meta, sep='\t', header=0)  
    if args.column not in meta_data.columns.tolist():
        sys.exit(f'{args.column} not found in meta file')
    df1 = meta_data[['accession', args.column]].copy()
    fasta_names = accessionMatch(args.fasta, meta_data)
    df2 = df1[df1['accession'].isin(fasta_names)].copy()
    species_dict = taxonDict()
    df2['tax_id'] = df2[args.column].map(species_dict)
    df3 = df2[['accession', 'tax_id']]
    df3["tax_id"] = df3["tax_id"].astype("Int64")
    print(df3)
    df3.to_csv(args.output, sep='\t', index=False, header=False)
    print('\nNaNs:')
    print(df3[df3['tax_id'].isna()])
  
    
def accessionMatch(fasta_file, meta):
    names = set()
    for rec in SeqIO.parse(open(fasta_file, 'r'), 'fasta'):
        name = rec.id
        names.add(name)
    fasta_names = sorted(list(names))
    meta_names = sorted(list(set(meta['accession'].to_list())))
    not_found = []
    for name in fasta_names:
        if name not in meta_names:
            not_found.append(name)
    if len(not_found) == 0:
        printInfo('Accessions match in fasta and meta data')
    else:
        print('[ERROR] Different accessions in meta data and fasta file')
        for item in not_found:
            print(item)
        sys.exit('Stopped')
    return fasta_names

def main():
    args = parse_args()
    args = checkArgs(args)
    printMessage('Assigning taxids to accessions')
    dbTaxonomy(args)
    print('Done')

if __name__ == "__main__":
   main()   

