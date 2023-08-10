'''Idenitfy the top K-mers (with K=15, say),
where the value of each K-mer is given by coverage
of the epimers (e-mers with e=9, typically)'''

import math
import argparse
from collections import Counter
from functools import cache
import random
from types import MappingProxyType as FrozenDict

import verbose as v
from bestsofar import BestSoFar
import sequtil
import episequtil

def getargs():
    '''get arguments from command line'''
    arg_parser = argparse.ArgumentParser(description=__doc__)
    paa = arg_parser.add_argument

    paa("-E",type=int,default=9,
        help="Epitope length")
    paa("-K",type=int,default=15,
        help="Big-k K-mer length")
    paa("-N",type=int,default=10,
        help="Compute top N kmers (including any kmers from --initseq)")
    paa("--inseq","-i",
        help="input sequence filename")
    paa("--initseq",
        help="Sequence file with initial epigraph sequences")
    paa("--rounds","-R",type=int,default=0,
        help="iterative rounds of randomized improvements")
    paa("--tmr",type=int,default=0,
        help="(too many repeats) remove sequences with more than TMR repeats")
    paa("--reverse",action="store_true",
        help="reverse sequences")
    paa("--reduce",action="store_true",
        help="reduce the pool of kmers by a factor of 1/4")
    paa("--seed",type=int,default=0,
        help="random number seed")
    paa("--output","-o",
        help="write output to file")
    paa('--verbose','-v',action="count",default=0,
        help="verbose")
    args = arg_parser.parse_args()
    return args

@cache
def get_mers_from_seq(k,seq,stride=1):
    '''
    from an input string, return a /set/ of distinct k-mers
    eg, k=3, seq=ABCDEFABCG -> {ABC,BCD,CDE,DEF,EFA,FAB,BCG}
    '''
    N = math.ceil(1 + (len(seq)-k)/stride)
    return set(seq[n*stride:n*stride+k] for n in range(N))

def get_mers_from_seqlist(k,seqlist,stride=1):
    '''
    from list of input strings, return a single set of distinct k-mers
    that appear in the strings.
    '''
    ## typical use: convert k-mers into e-mers
    emerset = set()
    for seq in seqlist:
        emersubset = get_mers_from_seq(k,seq,stride=stride)
        emerset.update( emersubset )
    return emerset

class MerManager:
    '''class for dealing with E-mers (epimers) and K-mers'''
    def __init__(self,K,E):
        assert K > E
        self.K = K
        self.E = E
        ## initialized separately (with set_... functions, below)
        self.kmers_pool = None
        self.epicount = None

    def set_epicount(self,epicount):
        self.epicount = FrozenDict(epicount) ### make an immutable copy

    def set_kmers_pool(self,kmers):
        self.kmers_pool = kmers

    def kmers_from_seq(self,seq):
        stride = self.K - self.E + 1
        return get_mers_from_seq(self.K,seq,stride=stride)

    def all_kmers_from_seq(self,seq):
        return get_mers_from_seq(self.K,seq,stride=1)

    def epimers_from_seq(self,seq):
        return get_mers_from_seq(self.E,seq,stride=1)

    def epimers_from_kmers(self,kmers):
        return get_mers_from_seqlist(self.E,kmers)

    def kmer_count(self,kmers_pool=None,epicount=None):
        '''
        score each kmer based on sum of score of epimers in kmer;
        return dict (actually Counter()) with score for each kmer
        '''
        epicount = epicount or self.epicount
        kmers_pool = kmers_pool or self.kmers_pool
        kmer_counter = Counter()
        for kmer in kmers_pool:
            kmer_counter[kmer] = sum(epicount[epi]
                                     for epi in self.epimers_from_seq(kmer))
        return kmer_counter

    def eval_kmers(self,kmers,*xtrakmers,epicount=None):
        '''return the score associated with a list (or set, etc) of kmers.
        note: can have multiple lists on the command line, they are all
        merged into one big list: eg, m_mgr.eval_kmers(kmers_init,kmers_new)
        '''
        epicount = epicount or self.epicount
        kmers = list(kmers) ## make a copy
        if xtrakmers:
            for xkmers in xtrakmers:
                kmers.extend(list(xkmers))
        return sum(epicount[epi] for epi in self.epimers_from_kmers(kmers))

## next fcns, they could be MerManager methods
## but for now are separate functions

def get_best_kmers(m_mgr,N,kmers_pool=None,
                   kmers_init=None,epi_counts=None):
    '''return list of next N kmers, greedily chosen to maximize epimer coverage'''

    ## we'll be setting values to zero, so only do that to the copy
    kmers_pool = kmers_pool or m_mgr.kmers_pool
    epi_counts = epi_counts or m_mgr.epicount
    epi_counts = Counter(epi_counts) ## make a (mutable) copy
    kmers_init=list(kmers_init) if kmers_init else []

    assert len(kmers_init) <= N

    tsum = m_mgr.eval_kmers(kmers_init)
    next_kmers=[]
    for n in v.vtqdm(range(N)):
        if n < len(kmers_init):
            top_kmer = kmers_init[n]
            tval = 0
        else:
            big_k_counter = m_mgr.kmer_count(epicount=epi_counts)
            top_kmer = max(big_k_counter,key=big_k_counter.get)
            tval = big_k_counter[top_kmer]

        tsum += tval
        v.vvprint("%3d %s %6d %7d" % (n+1,top_kmer,tval,tsum))
        for epi in m_mgr.epimers_from_seq(top_kmer):
            epi_counts[epi] = 0
        next_kmers.append(top_kmer)

    ## At this point, can we evaluate most valuable epimers we missed
    remaining_epi = sorted(epi_counts,key=epi_counts.get,reverse=True)
    v.vprint("Remaining epimers:")
    for n,epi in enumerate(remaining_epi):
        if n>10:
            break
        v.vprint("  ",epi,epi_counts[epi])

    return tsum,next_kmers

def restart_kmers(m_mgr,N,kmers_pool=None,fix=0.5):
    '''
    random initialize fraction (fix) of the kmers, and
    with those kmers fixed, optimize the rest;
    and then fix the optimized kmers, and
    re-optimize what were random initializers
    '''
    kmers_pool = kmers_pool or m_mgr.kmers_pool
    num_fixed = math.ceil(fix*N)
    kmers_init = random.sample(kmers_pool,num_fixed)
    val,kmers = get_best_kmers(m_mgr,N,
                               kmers_init=kmers_init)

    kmers_init = kmers[num_fixed:]
    val,kmers = get_best_kmers(m_mgr,N,
                              kmers_init=kmers_init)
    return val,kmers

def tweak_kmers(m_mgr,kmers_current,trials=10,fix=0.9):
    '''further optimize a list of K-mers'''
    N = len(kmers_current)
    score = m_mgr.eval_kmers(kmers_current)
    sofar = BestSoFar(score,kmers_current)
    for _ in range(trials):
        kmers = sofar.best
        ## fix fraction of the K-mers
        fix_kmers = random.sample(kmers,math.ceil(fix*N))
        val,kmers = get_best_kmers(m_mgr,N,kmers_init=fix_kmers)
        sofar.update(val,kmers)
    return val,kmers

def main(args):
    '''main: get top K-mers'''
    v.vprint(args)
    if args.seed:
        random.seed(args.seed)

    episequtil.init_epimer_length(args.E)

    seqs = sequtil.read_seqfile(args.inseq,rmdash=True,badchar='x')
    seqs = episequtil.seqgen_filter(seqs,reverse=args.reverse,
                                    tmr=args.tmr,rmx=False)

    mer_manager = MerManager(args.K,args.E)

    kmers = set()
    epimer_counter = Counter()
    for s in seqs:
        epimer_counter.update( mer_manager.epimers_from_seq(s.seq) )
        kmers.update( mer_manager.all_kmers_from_seq(s.seq) )

    mer_manager.set_epicount(epimer_counter)

    if args.reduce:
        ## reduce pool of kmers to those with higher scores
        counter = mer_manager.kmer_count(kmers_pool=kmers)
        top_kmers = sorted(counter,key=counter.get,reverse=True)
        kmers = set(top_kmers[:len(top_kmers)//4])
        v.print('kmers:',len(top_kmers),"reduced to",len(kmers))

    ## if there are initial sequences, then acquire initial k-mers from them
    kmers_init = set()
    if args.initseq:
        seqs = sequtil.read_seqfile(args.initseq,rmdash=True)
        seqs = episequtil.seqgen_filter(seqs,reverse=args.reverse)
        for s in seqs:
            v.vvprint('init:',s.name)
            kmers_init.update( mer_manager.kmers_from_seq(s.seq) )
        v.vprint('initialized with',len(kmers_init),'kmers')
    kmers.update(kmers_init)
    kmers=list(kmers)
    v.vprint('total kmers:',len(kmers))
    
    mer_manager.set_kmers_pool(kmers)

    ## All set up to get the rest of the kmers
    tsum,nkmers = get_best_kmers(mer_manager,args.N,
                                 kmers_init=kmers_init)

    v.vprint('next kmers:',len(nkmers),len(kmers_init))

    score = mer_manager.eval_kmers(nkmers)
    print(tsum,score)

    sofar = BestSoFar(score,nkmers)
    print("value:",sofar.value)

    ## let's try some random s..ff
    for _ in range(args.rounds):
        val,nkmers = restart_kmers(mer_manager,args.N)
        v.vprint('rnd_init:',val,'vs best:',sofar.value)

        val,nkmers = tweak_kmers(mer_manager,nkmers,trials=1,fix=0.5)
        print('new_kmer:',val,'vs best:',sofar.value)

        val,nkmers = tweak_kmers(mer_manager,nkmers,fix=0.9)
        print('nxt_kmer:',val,'vs best:',sofar.value)

        sofar.update(val,nkmers)
        nkmers = sofar.best
        val,nkmers = tweak_kmers(mer_manager,nkmers,fix=0.9)
        print('xxt_kmer:',val)

        sofar.update(val,nkmers)

    v.print("Score:",sofar.value)
    if args.output:
        ## should output be a .fasta? (it's a .seq)
        with open(args.output,'w') as f_out:
            for kmer in sofar.best:
                print(kmer,file=f_out)

if __name__ == '__main__':
    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
