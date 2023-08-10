''' epiun: EPIgraph UNaligned '''

# python epiun.py -E9 --inseq epi_inseqs_18June2019.tbl --nodelowcount 1 -v -M 2

## T==1 corresponds to a single trial, best and best-complement
## T>1, multiple trials, random and best-complement
##
## Mx==0: corresponds to "1+1", ie best and best-complement
## Mx>0: corresponds to "1+1+". ie iterative refinement
##       1. c[0]=best, then c[1]=best-complement, ..., c[M-1]=best-complement
##       2. Toss c[0]
##       3. Let c[0..M-2] <= c[1..M-1],
##       4. Compute c[M-1] = best complement to c[0..M-2]
##       5. Go to 2, up to Mx times
## Can make Mx >> 0, it will stop automatically when it converges

## Standard python imports
import re
import argparse
import random
import pickle
from collections import deque
import numpy as np # used for statistics of results
import networkx as nx

## Local not-epi-specific imports
import verbose as v
import sequtil

## Local imports
import config
import episequtil
import epiuutil as eu

import nxutil
Begin = nxutil.Begin
End   = nxutil.End

def getargs():
    '''parse the command line for options and arguments'''
    parser = argparse.ArgumentParser(description=__doc__)
    paa=parser.add_argument
    paa("-E",type=int,default=9,
        help="Epitope length")
    paa("--inseq","-i",
        help="input sequence filename")
    paa("--out",help="output file with vaccine sequences")
    paa("--outseqnames",
        help="comma-separated list of names of vaccine sequences")
    paa("-T",type=int,default=1,
        help="Number of trials")
    paa("--seed",type=int,default=0,
        help="random seed")
    paa("--Mxtra",type=int,default=0,
        help="extra iterations beyond M")

    paa("-M",type=int,default=1,
        help="Number of seqs in cocktail")

    paa("--weightsfile",
        help="File with weights associated to names")

    ### Not in principle needed, but might nudge a higher score
    paa("--tmr",type=int,default=0,
        help="(too many repeats) remove sequences with more than TMR repeats")
    paa("--reverse",action="store_true",
        help="reverse sequences")
    paa("--nodelowcount",type=int,default=0,
        help="delete nodes with count <= LOWCOUNT")
    paa("--edgelowcount",type=int,default=0,
        help="delete edges whose nodes have low counts")
    paa("--rmweaknodes",action="store_true",
        help="delete nodes whose self and all neighbors have count=1")

    ### Graph-specific
    paa("--dag",default="rss",
        help="Heuristic for dag-ifying graph: rss, sum, etc")
    paa("--writegraph",help="write graph object to python pickle file")
    paa("--usegraph",help="read graph object from python pickle file")

    ## Initialize
    paa("--fixinit",help="file with fixed initial vaccine sequences")
    paa("--hintinit",help="file with suggested initial vaccine sequences")

    paa("--verbose","-v",help="verbose",action="count",default=0)
    args = parser.parse_args()
    return args

def getinitseqs(filename,reverse=False):
    if filename:
        seqs = list(sequtil.read_seqfile(filename,rmdash=False))
        seqs = list(episequtil.seq_filter(seqs,reverse=reverse))
    else:
        seqs = []
    return seqs

def initcocktail(seqs):
    return [ [Begin] + eu.episeq_from_string(s.seq) + [End]
             for s in seqs
    ]

def graph_from_pkl(pklfilename,reverse=False):
    EPIMER = config.EPIMER
    with open(pklfilename,'rb') as G_infile:
        G = pickle.load(G_infile)
        G_topnodes = pickle.load(G_infile)
    v.vprint(f"Loaded graph/topnodes from {G_infile} with "
           f"{len(G)} = {len(G_topnodes)} nodes")
    ## Check this G is consistent with -E and --reverse
    G_EPIMER = G.graph.get('EPIMER',0) ## if not specified, don't trust it
    if G_EPIMER != EPIMER:
        raise RuntimeError(f"Graph [{pklfilename}] built with "
                           f"E={G_EPIMER} "
                           f"But -E{EPIMER} specified here")
    G_reverse = G.graph.get('reverse',False) ## if not specified, assume False
    if G_reverse != reverse:
        raise RuntimeError(f"Graph [{pklfilename}] built with "
                           f"reverse={G_reverse} "
                           f"But --reverse={reverse} specified here")

    return G,G_topnodes

def wacky_initialize(G,goodpath):
    ## Choose a node from middle of the good path (avoid Begin,End)
    midposn = random.randrange(1,len(goodpath)-1)
    newpath=[]
    if random.randrange(0,2) == 1:
        ## Good path from begin to mid, new path from mid to end
        newpath.extend(goodpath[:midposn])
        newpath.extend(nxutil.random_path_to_end(G,goodpath[midposn]))
    else:
        ## New path from begin to mid, then good path mid to end
        newpath.extend(nxutil.random_path_from_begin(G,goodpath[midposn]))
        newpath.extend(goodpath[midposn+1:])
    return newpath

def var_cocktail_init(t,G,
                      hintinit_seqs=None,
                      bestcocktail=None):

    var_cocktail = deque()

    if t==0 and hintinit_seqs is not None:
        ## On first trial, initialize with hint sequences
        for episeq in initcocktail(hintinit_seqs):
            var_cocktail.append(episeq)

    if t>0:
        ## On subsequent trials, initialize with random path
        if t % 2 == 1 or bestcocktail is None or len(bestcocktail)==0:
            ## On odd trials, just pick a random path
            var_cocktail.append(nxutil.random_path(G))
        else:
            ## On even trials, try something wacky
            goodpath = random.choice(bestcocktail)
            newpath = wacky_initialize(G,goodpath)
            var_cocktail.append(newpath)

    return var_cocktail

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
    config.EPIMER=EPIMER=args.E

    if args.M==1 and args.T>1:
        raise RuntimeError("If M=1, then T>1 not meaningful; use T=1")
    if args.seed:
        random.seed(args.seed)
        np.random.seed(args.seed)

    weights_byname = read_weights_file(args.weightsfile)

    ## GET GRAPH (from pickle file, or build one from sequence file)

    if args.usegraph:
        G,G_topnodes = graph_from_pkl(args.usegraph,reverse=args.reverse)
    else:
        seqs = sequtil.read_seqfile(args.inseq,rmdash=True,badchar='x')
        seqs = episequtil.seq_filter(seqs,reverse=args.reverse,
                                  tmr=args.tmr,rmx=False)
        G = eu.make_graph(seqs,dagify=args.dag,
                          nodelowcount=args.nodelowcount,
                          edgelowcount=args.edgelowcount,
                          rmweaknodes=args.rmweaknodes,
                          weights=weights_byname)
        G.graph['EPIMER']=EPIMER
        G.graph['reverse']=args.reverse
        v.vprint("Topological sorting of",len(G),"nodes")
        G_topnodes=nx.topological_sort(G) ## do this once, it won't change
        G_topnodes=list(G_topnodes) ## make a list, will use multiple times

    epi_total = G.graph['epiTotal']

    if args.writegraph:
        with open(args.writegraph, 'wb') as G_outfile:
            pickle.dump(G,         G_outfile, pickle.HIGHEST_PROTOCOL)
            pickle.dump(G_topnodes,G_outfile, pickle.HIGHEST_PROTOCOL)
        v.vprint("Dumped graph/topnodes to",G_outfile)

    if not nx.has_path(G,Begin,End):
        raise RuntimeError("Graph is broken!")

    ## OK, GRAPH HAS BEEN BUILT -- TIME TO FIND PATHS THROUGH IT

    ## get initial vaccine sequences/names from seq files
    fixinit_seqs  = list(getinitseqs(args.fixinit, reverse=args.reverse))
    hintinit_seqs = list(getinitseqs(args.hintinit,reverse=args.reverse))

    num_fixinit=len(fixinit_seqs)
    v.print('num_fixinit:',num_fixinit)
    if args.M <= num_fixinit:
        v.vprint(f"Too many --fixinit seqs: {num_fixinit}; "
               f"should be less than M={args.M}")

    num_hintinit=len(hintinit_seqs)
    if args.M-num_fixinit < num_hintinit:
        raise RuntimeError(f"Too many --hintinit sequences: {num_hintinit}; "
                           f"Should be no larger than "
                           f"M-num_fixinit={args.M}-{num_fixinit}"
                           f"={args.M-num_fixinit}")

    ## get names for the vaccine sequences
    ## fixed first, command-line-specifed second, assigned EG-# third
    vaccine_sequence_names = [ s.name for s in fixinit_seqs ]
    v.vvprint("Initial seq names:",vaccine_sequence_names)
    if args.outseqnames:
        vaccine_sequence_names.extend(args.outseqnames.split(','))
    for k in range(len(vaccine_sequence_names),args.M):
        vaccine_sequence_names.append("EG-%d"%k)
    v.vprint("Vaccine sequence names:",vaccine_sequence_names)
    assert len(vaccine_sequence_names) ==  args.M

    ## now that we have sequence names, we can write
    ## vaccine cocktails to the output file; here's the fcn for that
    def c_to_file(c):
        if len(c) > args.M:
            raise RuntimeError("Cocktail has more [%d] sequences, but M=%d"%
                               (len(c),args.M))
        episequtil.cocktail_to_file(c,args.out,
                                    seqnames=vaccine_sequence_names,
                                    reverse=args.reverse)

    ## and here's the fcn for computing cocktail coverage
    def c_coverage(c):
        return eu.graph_cocktail_coverage(G,c)
        #return eu.cocktail_coverage(c,epicount)


    ## Fix the first num_fixinit antigens:
    fix_cocktail = list( initcocktail(fixinit_seqs) )
    assert num_fixinit == len(fix_cocktail)
    if num_fixinit:
        score = c_coverage(fix_cocktail)
        print(f"fix score={score} ({score/epi_total:.6f}) ")


    bestcocktail = [] ## list of paths
    bestscore = -1
    scorelist = []

    ## Loop over multiple trials of cocktail finding
    for t in range(args.T):
        ## a cocktail has a fixed component and a variable component
        var_cocktail = var_cocktail_init(t,G,
                                         hintinit_seqs=hintinit_seqs,
                                         bestcocktail=bestcocktail)

        cocktail = fix_cocktail + list(var_cocktail)

        iterates_with_no_gain = -1
        for m in range(num_fixinit, args.M + args.Mxtra):

            if iterates_with_no_gain >= args.M-1:
                ## M iterations with no gain -- we're stuck!
                break

            ## If more than args.M in cocktail, get rid of extras
            ## eg, if fixinit + hintinit > args.M
            while len(var_cocktail) > args.M - num_fixinit:
                print("popleft! t=",t)
                var_cocktail.popleft()

            cocktail = fix_cocktail + list(var_cocktail)
            current_coverage = c_coverage(cocktail)

            ## Now get rid of the oldest path to make room for the new
            lastpath=None
            if len(cocktail) == args.M:
                lastpath = var_cocktail.popleft()
            cocktail = fix_cocktail + list(var_cocktail)

            newpath = eu.optimal_complementary_path(G,G_topnodes,cocktail)

            var_cocktail.append(newpath)
            cocktail = fix_cocktail + list(var_cocktail)
            new_coverage = c_coverage(cocktail)

            if new_coverage > current_coverage:
                ## if gain, then reset the no_gain counter
                iterates_with_no_gain = 0
            else:
                ## if there was no gain, update 'no_gain' counter
                ## then pop that newpath and put back the old lastpath
                iterates_with_no_gain += 1
                var_cocktail.pop()
                if lastpath:
                    var_cocktail.append(lastpath)
                cocktail = fix_cocktail + list(var_cocktail)

        score = c_coverage(cocktail)

        scorelist.append(score)
        print(f"t={t} score={score} ({score/epi_total:.6f}) "
              f"old best score= {bestscore} ({bestscore/epi_total:.6f})")
        if score > bestscore:
            bestscore=score
            bestcocktail=cocktail
            c_to_file(bestcocktail)


    print(f"best cocktail score: {bestscore} = {c_coverage(bestcocktail)}")

    if args.T > 1:
        scorelist = np.array(scorelist)
        print("Mean score: %9.2f" % np.mean(scorelist))
        print("Stdv score: %9.2f" % np.std(scorelist))
        print("Min  score: %6.0f" % np.min(scorelist))
        print("Med  score: %6.0f" % np.median(scorelist))
        print("Max  score: %6.0f" % np.max(scorelist))


    c_to_file(bestcocktail)

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
