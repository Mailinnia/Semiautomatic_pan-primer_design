'''
Primer check
'''

import pandas as pd
from os import path
import os
import sys
from Bio.Seq import Seq
import math
import itertools
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
from decimal import Decimal, ROUND_HALF_UP
import subprocess

KT_TAXONOMY_DB = '<USER-PATH>'

def roundHalfUp(n):
    return int(Decimal(n).quantize(Decimal("1"), rounding=ROUND_HALF_UP))

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

def targetCount(db):
    db_path = path.dirname(db)
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

################# Blast #################
def expandDegeneratePrimers(primer_seq):
    """Generate all possible non-degenerate sequences for a degenerate primer."""
    bases = [ambiguous_dna_values[base] for base in primer_seq]  # Get possible bases for each position
    all_variants = [''.join(p) for p in product(*bases)]  # Generate all combinations
    return all_variants

def createBlastQuery(df1, query_out):
    print('Creating blast query')
    records = []
    for idx, row in df1.iterrows():
        primer_name = row['primer_name']
        primer_seq = row['seq']
        expanded_primers = expandDegeneratePrimers(primer_seq)
        for i, seq in enumerate(expanded_primers):
            name = f"{primer_name}_var{i+1}"
            rec = SeqRecord(Seq(seq), id=name, name=name, description=f"{primer_name} variant {i+1}")
            records.append(rec)
    #        
    with open(query_out, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    print(f'Created {len(records)} permutations')

def blastPrimers(db,  blast_out, n_seqs, blast_query, n_threads):
    cmd = ['blastn','-num_threads', n_threads, '-query', blast_query, '-db', db, '-out', blast_out, '-task', 'blastn-short', 
            '-outfmt', '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send sstrand qseq sseq staxid', '-evalue', '5000', 
            '-reward', '1', '-penalty', '-1', '-gapopen', '2', '-gapextend', '1', '-max_target_seqs', str(n_seqs)]
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
    process = subprocess.run(cmd)
    if process.returncode !=0:
        sys.exit(f'{cmd} FAILED')

################# Taxon dict #################
def taxonDict():
    print('Preparing taxonomies')
    taxon_db = pd.read_csv(KT_TAXONOMY_DB, sep = '\t', header=0, names=['taxid', 'node', 'parent', 'level', 'name'])
    taxon_dict = taxon_db[['taxid', 'level', 'parent', 'name']].set_index('taxid').to_dict('index')
    taxon_dict[0]={'level': 'no rank', 'parent': 1, 'name': 'unclassified'}
    print('  Finished preparations')
    return taxon_dict

def traceTaxa(taxon_dict, key):
    if key not in taxon_dict:
        if key == 0:
            return [('no rank', 'unclassified', 0, 1)]
        else:
            return []
    taxon = taxon_dict[key]
    parent = taxon['parent']
    if parent == 1:
        return [(taxon['level'], taxon['name'], key,taxon['parent'])]
    return [(taxon['level'], taxon['name'], key,taxon['parent'])] + traceTaxa(taxon_dict, parent)

def getSpecies(taxon_dict, key):
    trace = traceTaxa(taxon_dict, key)
    for lvl, name, _, _ in trace:
        if lvl =='species':
            return name
        

################# On-targets #################
def matchAmbig(base1, base2):
    return bool(set(ambiguous_dna_values.get(base1, base1)) & set(ambiguous_dna_values.get(base2, base2)))

def match3prime(row):
    query_end = row["query_seq"][-3:]
    ref_end = row["ref_seq"][-3:]
    return all(matchAmbig(q, r) for q, r in zip(query_end, ref_end))

def countAmbigs(row):
    query_seq = row["query_seq"]
    ref_seq = row["ref_seq"]
    ambig = ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']
    count = 0
    for q, r in zip(query_seq, ref_seq):
        if q == r:
            continue
        if q in ambig or r in ambig:
            count += 1
    return count

def filterTargetBlast(blast_out, filter_primers, target_mismatch):
    """
    create a BLAST hit database for each primer. filter for mismatches.
    returns a prefiltered pandas df
    """
    print(f'Target mismatch set to: {target_mismatch}')
    naming_dict = filter_primers.set_index('primer_name')[['amp_name', 'oligo_type','amp_len']].to_dict() 
    #
    columns = ["query", "ref", "query_len", "aln_len", "mismatch", "gaps", 'query_start', 'query_end', "ref_start", "ref_end","strand" ,'query_seq', 
                'ref_seq', 'tax_id']
    df1 = pd.read_table(blast_out, names=columns)
    # Extract the last accession number from the 'ref' column
    # df1['ref'] = df1['ref'].str.extract(r'(\w+\|\w+\.\d+\|)$')
    # Remove the leading database prefix (e.g., "gb|", "dbj|", "emb|") to keep only the accession number
    # df1['ref'] = df1['ref'].str.extract(r'(\w+\.\d+)')
    accessions = df1['ref'].str.extract(r'(?:\w+\|)?(\w+\.\d+)(?:\|)?$')[0]
    df1['ref'] = accessions.combine_first(df1['ref'])
    df1['query'] = df1["query"].str.replace(r"_var\d+$", "", regex=True)
    print(df1)
    df1['amp_name'] = df1['query'].map(naming_dict['amp_name'])
    df1['oligo_type']=df1['query'].map(naming_dict['oligo_type'])
    df1['amp_len']=df1['query'].map(naming_dict['amp_len'])
    df2 = df1[(df1['query_end'] == df1['query_len']) & # must match in 3' end of primer
                ((df1["query_len"] - df1["aln_len"]) <= round(df1["query_len"] * 0.25))].copy() # must align to majority of primer
    df3 = df2[((df2["oligo_type"] == "left") & (df2["strand"] == "plus")) |  # Keep left primers that are on plus strand
            ((df2["oligo_type"] == "right") & (df2["strand"] == "minus")) |  # Keep right primers that are on minus strand
            (df2["oligo_type"] == "probe")  # Keep probes regardless of strand 
                ].copy()
    df3["ambiguity_count"] = df3.apply(countAmbigs, axis=1)
    df3["adjusted_mismatch"] = df3["mismatch"] - df3["ambiguity_count"] #ambigs count as mismatches. So we adjust for those
    df4 = df3[df3.apply(match3prime, axis=1)] # no ambiguities in the first 3 pos of 3'
    df5 = df4[df4['adjusted_mismatch'] <= target_mismatch].copy() # remove rows with mismatches above target_mismatch
    df5 = df5.sort_values(['query', 'ref', 'adjusted_mismatch'])
    df6 = df5.drop_duplicates(subset=['query', 'ref', 'ref_start', 'ref_end'], keep='first').reset_index(drop=True)
    tax_ids = list(set(df6['tax_id'].fillna('NaN').to_list()))
    if len(tax_ids) == 1 and tax_ids[0]=='NaN':
        df6['species'] = 'n.d.'
    else:
        taxon_dict = taxonDict()
        df6['species'] = df6['tax_id'].map(lambda taxid: getSpecies(taxon_dict, taxid))

    return df6
    # removes tabular output
    # os.remove(blast_out)
    
def targetAmpDict(blast_df, ref_result, ref_set):
    amp_dict = {'success':{}, 'fails':{}, 'fail_species':{}}
    species_dict = blast_df.set_index('ref')['species'].to_dict()
    for name in ref_result.keys():
        name_ref = ref_result[name]
        fail_acc = []
        fail_species = []
        success = []
        for key, value in name_ref.items():
            if value is False:
                fail_acc.append(key)
                species = species_dict[key]
                if species is None:
                    print(f'FOUND NONE:{key}')
                fail_species.append(species)  
            if value is True:
                success.append(key)
        refs = set(name_ref.keys())
        refs_missing = list(ref_set-refs)
        fail_acc = fail_acc+refs_missing
        fail_str = '; '.join(fail_acc)
        if fail_str == '':
            fail_str = 'None'
        fail_spec_str = '; '.join(sorted(set(x for x in fail_species if x is not None)))
        amp_dict['success'][name]=len(success)
        amp_dict['fails'][name]=fail_str
        amp_dict['fail_species'][name]=fail_spec_str
    #
    return amp_dict

def targetRange(amp_len):
    ranges = {i: i // 10 for i in range(1000, 20001, 500)}
    ranges[300]=30
    for key in sorted(list(ranges.keys())):
        if amp_len <= key:
            return ranges[key]


def checkTargets(filter_primers, blast_df, amp):
    filter_amp = filter_primers[filter_primers['amp_name']==amp]
    oligo_types = set(filter_amp['oligo_type'].to_list())
    #
    left_name = filter_amp[filter_amp['oligo_type']=='left']['primer_name'].to_list()[0]
    right_name = filter_amp[filter_amp['oligo_type']=='right']['primer_name'].to_list()[0]
    ref_result = {left_name:{}, right_name:{}} #initiate dict
    #
    probe_present = False
    if 'probe' in oligo_types:
        probe_present = True
        probe_names = set(filter_amp[filter_amp['oligo_type']=='probe']['primer_name'].to_list())
        for probe in probe_names:
            ref_result[probe] = {} #initiate dict with probes
    #
    ref_set = set(blast_df['ref'].to_list())
    df1 = blast_df[blast_df['amp_name']==amp]
    amp_len = df1['amp_len'].to_list()[0]
    #
    for ref in set(df1['ref']):
        check_probes = probe_present
        df2 = df1[df1['ref']== ref]
        left_primers = df2[df2['oligo_type']=='left']['ref_start'].to_list()
        right_primers = df2[df2['oligo_type']=='right']['ref_start'].to_list()
        if len(left_primers) ==0 or len(right_primers) == 0:
            ref_result[left_name][ref] = False
            ref_result[right_name][ref] = False
            if check_probes: 
                for probe in probe_names:
                    ref_result[probe][ref] = False
            continue #skip this ref
        combinations = [(l, r, r - l) for l, r in itertools.product(left_primers, right_primers)]
        #
        if probe_present:
            probe_starts = df2[df2['oligo_type']=='probe']['ref_start'].to_list()
            probe_ends = df2[df2['oligo_type']=='probe']['ref_end'].to_list()
            blast_probes = df2[df2['oligo_type']=='probe']['query'].to_list()
            if len(probe_starts) ==0:
                for probe in probe_names:
                    ref_result[probe][ref] = False
                check_probes = False #no need to check placements of probes
        #
        amp_ok = False
        target_diff = targetRange(amp_len)
        for left_start, right_end, combi_len in combinations:
            if combi_len >= amp_len-target_diff and combi_len <= amp_len+target_diff:
                # print(f'{left_start} and {right_end} within size limit')
                if check_probes:
                    for probe_start, probe_end, blast_probe in zip(probe_starts, probe_ends, blast_probes):
                        if probe_start >= left_start and probe_end <= right_end:
                            ref_result[blast_probe][ref] = True
                amp_ok = True
                break
        #
        ref_result[left_name][ref] = amp_ok
        ref_result[right_name][ref] = amp_ok
        if probe_present:
            for probe in probe_names:
                if ref not in ref_result[probe]:
                    ref_result[probe][ref]=False
    #
    amp_dict = targetAmpDict(blast_df, ref_result, ref_set)
    #
    return amp_dict     

def predictCatchTargets(filter_primers, blast_df, n_threads):
    '''
    Main function to predict whether all targets are caught by the designed primers
    '''
    amplicons = sorted(set(filter_primers['amp_name'].to_list()))
    with multiprocessing.Pool(processes=int(n_threads)) as pool:
        results = pool.starmap(
            checkTargets,
            [(filter_primers, blast_df, amp) for amp in amplicons]
        )
    target_dict = {'success': {}, 'fails': {}, 'fail_species':{}}
    for result in results:
        if result is not None:  # Ensure valid results
            for key in list(result.keys()):
                if key in result:
                    target_dict[key].update(result[key])  # 
    #   
    return target_dict

def searchTargets(filter_primers, n_seqs, blast_query, n_threads, db, blast_out, allowed_mismatch):
    print('\n[INFO] Searching for on-targets in provided blast-db')
    blastPrimers(db, blast_out, n_seqs, blast_query, n_threads)
    blast_df = filterTargetBlast(blast_out, filter_primers, allowed_mismatch)
    target_dict = predictCatchTargets(filter_primers, blast_df, n_threads)
    return target_dict


################# Off-targets #################
def count3primeMismatches(row, n=3):
    # Extract the last n bases from the query primer sequence
    # Extract the last n bases from the BLAST-aligned reference sequence
    # Compare the query and reference base-by-base: zip(query_end, ref_end) creates pairs: [('T','C'), ('C','C'), ('G','G')]
    # matchAmbig(q, r) returns: True  -> bases are compatible, False -> mismatch
    # "not matchAmbig(q, r)" converts: True  -> False  (0 mismatches), False -> True   (1 mismatch)
    # sum(...) counts how many mismatches occur
    query_end = row["query_seq"][-n:]
    ref_end = row["ref_seq"][-n:]

    return sum(
        not matchAmbig(q, r)
        for q, r in zip(query_end, ref_end)
    )

def filterOffTargetBlast(blast_out, filter_primers, allowance_3_prime):
    """
    create a BLAST hit database for each primer. filter for mismatches.
    returns a prefiltered pandas df
    """
    taxon_dict = taxonDict()
    naming_dict = filter_primers.set_index('primer_name')[['amp_name', 'oligo_type','amp_len']].to_dict() 
    #
    columns = ["query", "ref", "query_len", "aln_len", "mismatch", "gaps", 'query_start', 'query_end', "ref_start", "ref_end","strand" ,'query_seq', 
               'ref_seq', 'tax_id']
    df1 = pd.read_table(blast_out, names=columns)
    # # Extract the last accession number from the 'ref' column
    # df1['ref'] = df1['ref'].str.extract(r'(\w+\|\w+\.\d+\|)$')
    # # Remove the leading database prefix (e.g., "gb|", "dbj|", "emb|") to keep only the accession number
    # df1['ref'] = df1['ref'].str.extract(r'(\w+\.\d+)')    accessions = df1['ref'].str.extract(r'(?:\w+\|)?(\w+\.\d+)(?:\|)?$')[0]
    accessions = df1['ref'].str.extract(r'(?:\w+\|)?(\w+\.\d+)(?:\|)?$')[0]
    df1['ref'] = accessions.combine_first(df1['ref'])
    df1['query'] = df1["query"].str.replace(r"_var\d+$", "", regex=True)
    df1['amp_name'] = df1['query'].map(naming_dict['amp_name'])
    df1['oligo_type']=df1['query'].map(naming_dict['oligo_type'])
    df1['amp_len']=df1['query'].map(naming_dict['amp_len'])
    df2 = df1[(df1['query_end'] == df1['query_len']) & # must match in 3' end of primer
                ((df1["query_len"] - df1["aln_len"]) <= round(df1["query_len"] * 0.50))].copy() # must align to at least half of the primer
    df2["ambiguity_count"] = df2.apply(countAmbigs, axis=1)
    df2["adjusted_mismatch"] = df2["mismatch"] - df2["ambiguity_count"] #ambigs count as mismatches. So we adjust for those
    df2["mismatch_3prime"] = df2.apply(count3primeMismatches, axis=1) # count mismatches in the 3 bases of the 3' end
    df3 = df2[df2["mismatch_3prime"] <= allowance_3_prime].copy()  # allow up to n mismatch in the 3'-end, as this is off-target
    # df3 = df2[df2.apply(match3prime, axis=1)].copy() # no mismatches in the first 3 pos of 3'
    df3['allowed_mismatch'] = (df1["aln_len"]-3)/5 # allow mismatch 1 per 5 bp after the 3 first bases in 3' end
    df3['allowed_mismatch'] = df3['allowed_mismatch'].apply(lambda x: roundHalfUp(x)) # round half values up
    df4 = df3[df3['adjusted_mismatch'] <= df3['allowed_mismatch'] ].copy() # remove rows with mismatches above allowed mismatch
    df4 = df4.sort_values(['query', 'ref', 'adjusted_mismatch'])
    df5 = df4.drop_duplicates(subset=['query', 'ref', 'ref_start', 'ref_end'], keep='first').reset_index(drop=True)

    df5['species'] = df5['tax_id'].map(lambda taxid: getSpecies(taxon_dict, taxid))
    # removes tabular output
    # os.remove(blast_out)
    return df5

def checkOffTargets(filter_primers, blast_df, amp):
    filter_amp = filter_primers[filter_primers['amp_name']==amp]
    oligo_types = set(filter_amp['oligo_type'].to_list())
    #
    left_name = filter_amp[filter_amp['oligo_type']=='left']['primer_name'].to_list()[0]
    right_name = filter_amp[filter_amp['oligo_type']=='right']['primer_name'].to_list()[0]
    ref_result = {left_name:{}, right_name:{}} #initiate dict
    #
    probe_present = False
    if 'probe' in oligo_types:
        probe_present = True
        probe_names = set(filter_amp[filter_amp['oligo_type']=='probe']['primer_name'].to_list())
        for probe in probe_names:
            ref_result[probe] = {} #initiate dict with probes
    #
    df1 = blast_df[blast_df['amp_name']==amp]
    amp_len = df1['amp_len'].to_list()[0]
    min_len = math.floor(amp_len*0.5)
    max_len = max(1000, math.ceil(amp_len*1.5))
    #
    for ref in set(df1['ref']):
        check_probes = probe_present
        df2 = df1[df1['ref']== ref]
        primers = df2[df2['oligo_type']!='probe']['ref_start'].to_list()
        combinations = list(itertools.combinations(primers, 2))
        #
        if probe_present:
            probe_starts = df2[df2['oligo_type']=='probe']['ref_start'].to_list()
            probe_ends = df2[df2['oligo_type']=='probe']['ref_end'].to_list()
            blast_probes = df2[df2['oligo_type']=='probe']['query'].to_list()
            if len(probe_starts) ==0:
                for probe in probe_names:
                    ref_result[probe][ref] = False
                check_probes = False #no need to check placements of probes
        #
        amp_ok = False
        for a, b in combinations:
            combi_len = abs(a-b)
            if combi_len >= min_len and combi_len <= max_len:
                # print(f'{left_start} and {right_end} within size limit')
                if check_probes:
                    for probe_start, probe_end, blast_probe in zip(probe_starts, probe_ends, blast_probes):
                        if probe_start >= min([a,b]) and probe_end <= max([a,b]):
                            ref_result[blast_probe][ref] = True
                amp_ok = True
                break
        #
        ref_result[left_name][ref] = amp_ok
        ref_result[right_name][ref] = amp_ok
        if probe_present:
            for probe in probe_names:
                if ref not in ref_result[probe]:
                    ref_result[probe][ref]=False
        #
    # amp_dict = dict()
    amp_dict = {'off_targets':{}, 'off_target_species':{}, 'off_targets_no':{}}
    species_dict = blast_df.set_index('ref')['species'].to_dict()
    for name in ref_result.keys():
        name_ref = ref_result[name]
        caught_acc = []
        caught_species = []
        for key, value in name_ref.items():
            if value is True:
                caught_acc.append(key)
                caught_species.append(species_dict[key])   
        # caughts = [key for key, value in name_ref.items() if value is True]
        caught_no = len(caught_acc)
        caught_str = '; '.join(caught_acc)
        caught_spec_str = '; '.join(sorted(set('None' if spec is None else spec for spec in caught_species)))

        if caught_spec_str == '':
            caught_spec_str = 'None'
        amp_dict['off_targets'][name]=caught_str
        amp_dict['off_target_species'][name]=caught_spec_str
        amp_dict['off_targets_no'][name]=caught_no
    return amp_dict     

def predictNonSpecific(filter_primers, blast_df, n_threads):
    '''
    Main function to predict whether all targets are caught by the designed primers
    '''
    amplicons = sorted(set(filter_primers['amp_name'].to_list()))
    with multiprocessing.Pool(processes=int(n_threads)) as pool:
        results = pool.starmap(
            checkOffTargets,
            [(filter_primers, blast_df, amp) for amp in amplicons]
        )
    target_dict = {'off_targets_no':{}, 'off_targets':{}, 'off_target_species':{}}
    for result in results:
        if result is not None:  # Ensure valid results
            for key in list(target_dict.keys()):
                if key in result:
                    target_dict[key].update(result[key])  # 
    #
    return target_dict

def searchOffTargets(filter_primers, n_seqs, allowance_3_prime, blast_query, n_threads, db, blast_out, exclude=None):
    print('\n[INFO] Searching for off-targets in provided blast-db')
    blastPrimers(db, blast_out, n_seqs, blast_query, n_threads)
    blast_df = filterOffTargetBlast(blast_out, filter_primers, allowance_3_prime)
    print(blast_df)
    offtarget_dict = predictNonSpecific(filter_primers, blast_df, n_threads)
    return offtarget_dict
