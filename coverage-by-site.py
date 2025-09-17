'''make plot of coverage by site of epigraphs '''

import re
import random
import argparse
from collections import Counter,defaultdict
import itertools as it
import matplotlib.pyplot as plt
import numpy as np

import verbose as v
import sequtil

import config ## defines globals; specifically EPIMER
import epiuutil as eu

DASH='-'
BADCHAR='x'

def getargs():
    '''get command line arguments'''
    parser = argparse.ArgumentParser()
    paa = parser.add_argument
    paa("--inseq","-i",
        help="input sequence filename")
    paa("--epigraphinput","-e",
        help="input epigraph filename")
    paa("-E",type=int,default=9,
        help="length of epitope")
    paa("-M",type=int,default=0,
        help="maximum number of epigraphs to consider")
    paa("--refpos",action="store_true",
        help="Use reference position (assume first seq is ref)")
    paa("--saveplot",
        help="write plot to file called SAVEPLOT")
    paa("--log",action="store_true",
        help="make plot semilogy")
    paa("--xrange",nargs=2,type=int,
        help="in plots, x-axis over this range")
    paa("--weightsfile","-w",
        help="tsv file assigns weights to sequences names")
    paa("-N",type=int,default=0,
        help="Only read N input sequences (for testing)")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = parser.parse_args()
    return args

def mkepitopecounter(seqstrings,dash=DASH, weights=None):
    '''From an aligned list of sequence strings
       make a count of epitopes at each position
       epicount[n][epi] is number of time epi is in position n

       in general, we worry about dashes
       For example (EPIMER=3)
                        long_epi   -> short_epi -> epi
       ...A-B---CDEF -> A-B---CDEF -> ABCDEF    -> ABC
       ...-B---CDDEF -> B---CDDEF  -> BCDDEF    -> -BC [*]
       ...A-B------- -> A-B------- -> AB        -> AB  [**]
       [*] this will be a "placeholder" epitope (count=0)
       [**] len != EPIMER, so not added to epilist
    '''

    EPIMER = config.EPIMER

    if not weights:
        weights = [1.0 for _ in seqstrings]

    epicount = defaultdict(Counter)
    if not dash:
        ## then it's easy!
        for seqstr,weight in zip(seqstrings,weights):
            for n in range(len(seqstr)-EPIMER+1):
                epi = seqstr[n:n+EPIMER]
                epicount[n][epi] += weight
        return epicount

    assert len(dash) == 1
    re_dash = re.compile(dash)
    for seqstr,weight in zip(seqstrings,weights):
        nd_seqstr = re_dash.sub('',seqstr)
        nd_ndx = 0
        nd_ndx_max = len(nd_seqstr)-EPIMER
        for n in range(len(seqstr)-EPIMER+1):
            if seqstr[n] == dash:
                continue
                #epi = dash + nd_seqstr[nd_ndx:nd_ndx+EPIMER-1]
            if nd_ndx > nd_ndx_max:
                break
            epi = nd_seqstr[nd_ndx:nd_ndx+EPIMER]
            nd_ndx += 1
            if (BADCHAR not in epi):
                epicount[n][epi] += weight

    return epicount

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
    '''main coverage-by-site'''
    v.vprint(args)

    config.EPIMER = args.E

    episeqs = sequtil.read_seqfile(args.epigraphinput,
                                   rmdash=True,badchar=BADCHAR)
    episeqs = list(episeqs)
    v.vprint(f"{len(episeqs)} epigraph sequences")
    if args.M:
        episeqs = episeqs[:args.M]

    allepitopes = []
    epitopes = set()
    for s in episeqs:
        epitopes.update( eu.episeq_from_string(s.seq) )
        allepitopes.append(epitopes.copy())
    v.vprint(f"{len(epitopes)} epitopes")

    seqs = sequtil.read_seqfile(args.inseq,
                                rmdash=False,badchar=BADCHAR)
    seqs = list(seqs)

    site_ndx=None
    if args.refpos:
        firstseq = seqs[0].seq
        site_ndx = it.accumulate((int(seqchar!=DASH)
                                  for seqchar in firstseq),initial=0)
        site_ndx = list(site_ndx)

    v.vprint(f"{len(seqs)} sequences read")
    if args.N:
        seqs = random.choices(seqs,k=args.N)
        v.vprint(f"{len(seqs)} sequences sampled")
    num_seqs = len(seqs)

    weights_list = None
    if args.weightsfile:
        weights_dict = read_weights_file(args.weightsfile)
        weights_list = [weights_dict[s.name] for s in seqs]
        num_seqs = sum(weights_list) ## effective number of seqs

    epicounts = mkepitopecounter([s.seq for s in seqs],
                                 weights=weights_list)


    if site_ndx:
        ## re-define epicount to include all epi's at a given reference position
        ## (some will be at different alignment positions; so add up all the
        ## values at all alignment positions corresponding to a reference position)
        amalg_epicounts = defaultdict(Counter) ## count epitopes at amalgamated ref site nn
        for n in sorted(epicounts):
            nn = site_ndx[n]
            for epi in epicounts[n]:
                amalg_epicounts[nn][epi] += epicounts[n][epi]
        epicounts = amalg_epicounts


    xyk=[[] for _ in allepitopes]

    ## loop over index of position in aligned sequence
    for n in sorted(epicounts):
        ## Loop over epigraph number; 'epitopes' is cumulative set covered
        tot = sum( epicounts[n].values() )
        for k,epitopes in enumerate(allepitopes):
            mis = sum( epicounts[n][epi]
                       for epi in epicounts[n] if epi not in epitopes )
            xyk[k].append((n,100*mis/num_seqs))


    for k,xy in enumerate(xyk):
        print(f"EPI-{k+1} ",end="")
        ## Overall:
        avemis = np.mean([mis for n,mis in xy])
        print(f"{avemis:7.4f} percent missed (approximate!)")

    ## make plot
    plot = plt.semilogy if args.log else plt.plot
    for k,xy in enumerate(xyk):
        x,y = zip(*xy)
        plot(x,y,label=f"EPI-{k+1}",lw=(k+1)/2)

    plt.ylabel("Missed epitopes (percent)")
    refalign = "Reference" if args.refpos else "Alignment"
    plt.xlabel(refalign + " Position")
    if args.xrange is not None:
        plt.xlim(args.xrange)
    if args.log:
        plt.ylim([1.0e-3,1.0e2])
    else:
        plt.ylim([0,75]) ## fixme -- hardcoded

    plt.legend(bbox_to_anchor=(1.02, 1),
               #frameon=False,
               #handletextpad=0,
               labelspacing=0.45,
               loc='upper left', borderaxespad=0.,
               prop={'family' : 'monospace'})

    plt.tight_layout()
    if args.saveplot:
        plt.savefig(args.saveplot)
    else:
        plt.show()

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
