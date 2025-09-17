'''compute coverage of epigraphs '''

import re
import random
import argparse
from collections import Counter,defaultdict
import itertools as it
import numpy as np

## Local not-epi-specific imports
import verbose as v
import sequtil

## Local imports
import config ## defines globals; specifically EPIMER
import episequtil
import epiuutil as eu

DASH='-'
BADCHAR='x'

def getargs():
    '''get command line arguments'''
    parser = argparse.ArgumentParser(description=__doc__)
    paa = parser.add_argument
    paa("-E",type=int,default=9,
        help="length of epitope")
    paa("--inseq","-i",
        help="input sequence filename")
    paa("--epigraphinput","-e",
        help="input epigraph filename")
    paa("--weightsfile","-w",
        help="tsv file assigns weights to sequence names")
    paa("--epicountfile",
        help="output tsv file with epitopes and total weight, for each epigraph")
    paa("-N",type=int,default=0,
        help="Only read N input sequences (for testing)")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = parser.parse_args()
    return args

def read_weights_file(file):
    '''
    read the weights file, make a weights dict
    each weight is weights[seq_name]
    '''
    catalog=dict()
    if not file:
        return catalog
    
    with open(file) as fin:
        for line in fin:
            line = re.sub('#.*','',line)
            line = line.strip()
            if not line:
                continue
            value,key = line.split()
            catalog[key]=float(value)
    return catalog

def main(args):
    '''main coverage'''
    v.vprint(args)

    config.EPIMER = args.E

    episeqs = sequtil.read_seqfile(args.epigraphinput,
                                   rmdash=True,badchar=BADCHAR)
    episeqs = list(episeqs)
    v.vprint(f"{len(episeqs)} epigraph sequences")
    per_episeq_epitopes = [] ## list of sets: k'th item is set of epitopes in k'th epigraph sequence
    epitopes_sofar = set()
    for s in episeqs:
        epitopes_in_episeq = set( eu.episeq_from_string(s.seq) )
        per_episeq_epitopes.append( epitopes_in_episeq )
        epitopes_sofar.update( epitopes_in_episeq )
        v.vprint(f'{len(epitopes_in_episeq)} epitopes in episeq; {len(epitopes_sofar)} epitopes so far')
    v.vprint(f"{len(epitopes_sofar)} epitopes total in all episeqs")
    
    seqs = sequtil.read_seqfile(args.inseq,
                                rmdash=False,badchar=BADCHAR)
    seqs = list(seqs)

    weights_list = [1]*len(seqs)
    if args.weightsfile:
        weights_dict = read_weights_file(args.weightsfile)
        weights_list = [weights_dict[s.name] for s in seqs]
        num_seqs = sum(weights_list) ## effective number of seqs
        v.vprint("Total number of sequences:",num_seqs)
        
    sum_mis = [0]*len(per_episeq_epitopes)
    sum_wt = 0
    cnt = [Counter() for _ in per_episeq_epitopes]
    for s,wt in zip(seqs,weights_list):
        sum_wt += wt
        epitopes_in_seq = set( eu.episeq_from_string(s.seq) )
        epi_missed = epitopes_in_seq.copy() ## start out missing them all
        for k,epitopes_in_episeq in enumerate(per_episeq_epitopes):
            epi_missed = epi_missed - epitopes_in_episeq
            mis = wt * len(epi_missed)/len(epitopes_in_seq)
            sum_mis[k] += mis
            if args.epicountfile:
                ## only need this info if writing the epicount files
                epi_covered = epitopes_in_seq & epitopes_in_episeq
                for epi in epi_covered:
                    cnt[k][epi] += wt

    for k,smis in enumerate(sum_mis):
        print(f'{k+1:2d} {100*smis/sum_wt:7.3f} percent missed')

    if args.epicountfile:
        for k,epitopes in enumerate(per_episeq_epitopes):
            with open(f'{args.epicountfile}-{k+1}.tsv','w') as fout:
                for epi in epitopes:
                    fout.write("%s\t%d\t%f\n" % (epi,int(cnt[k][epi]),cnt[k][epi]/sum_wt))
            
        
if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
