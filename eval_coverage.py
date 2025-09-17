'''Evaluate Coverage of a set of centroids on an input sequence list'''

import re
import random
import argparse
from collections import Counter
import numpy as np   ## used for statistics of bootstrap results

import verbose as v
import sequtil

import config ## defines globals; specifically EPIMER
import epiuutil as eu
import episequtil as epsu
import offby_count as oc

DASH='-'
BADCHAR='x'

def getargs():
    parser = argparse.ArgumentParser(description=__doc__)
    paa = parser.add_argument
    paa("-E",type=int,default=9,
        help="Epitope length")
    paa("--inseq","-i",
        help="input target population sequence filename")
    paa("--epigraphinput","-e",
        help="input epigraph (centroid/antigen/vaccine) filename")
    paa("-M",type=int,default=0,
        help="Number of antigens to read from vaccine file")
    paa("-B",type=int,default=0,
        help="Number of bootstraps to evaluate error")
    paa("--usebadepi",action="store_true",
        help="Put bad epitopes in denominator?")
    paa("--usecounts",action="store_true",
        help="Use counts (not fractions) for epitope coverage")
    paa("--rare",action="store_true",
        help="Include stats on rarest epitope")
    paa("--offby",type=int,default=0,
        help="Evaluate off-by-one(or two...) coverage")
    paa("--tmr",type=int,default=0,
        help="(too many repeats) remove strings with more than TMR repeats")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = parser.parse_args()
    return args


def epitope_tally_total(seqs,bad_epitopes_in_denominator=False):
    seq_epicount = eu.epitope_tally(seqs,rmx=False)
    if bad_epitopes_in_denominator:
        total = sum( seq_epicount.values() )
        epsu.remove_bad_epitopes(seq_epicount)
    else:
        epsu.remove_bad_epitopes(seq_epicount)
        total = sum( seq_epicount.values() )
    return total,seq_epicount

def main(args):
    """main eval_coverage"""
    
    config.EPIMER=args.E

    seqs = sequtil.read_seqfile(args.inseq,rmdash=True,badchar=BADCHAR)
    seqs = list(seqs)
    seqs = eu.seq_filter_repeats(seqs,args.tmr)

    cent = sequtil.read_seqfile(args.epigraphinput,rmdash=True)
    cent = list(cent)
    if args.M > len(cent):
        raise RuntimeError(f"fewer than C={arg.m} sequences"
                           f"in file: {args.epigraphinput}")
    if args.M <= 0:
        args.M = len(cent)
    if args.M < len(cent):
        cent = cent[:args.M]

    total,seq_epicount = epitope_tally_total(seqs,args.usebadepi)

    cen_epicount = eu.epitope_tally(cent)

    rogue_epi = set(cen_epicount.keys()) - set(seq_epicount.keys())
    for e in rogue_epi:
        seq_epicount[e]=0

    obarray = oc.mkobarray(seq_epicount.keys(),
                           cen_epicount.keys(),args.offby)

    stddevstr = " (std.dev.)" if args.B else ""
    header = "  m    exact"+stddevstr
    for d in range(1,args.offby+1):
        header = header + f" offby<={d}" + stddevstr
    if args.rare:
        header = header + " rarest rare-count"
    print(header)

    fmt   = " %8d"     if args.usecounts else " %8.6f"
    b_fmt = " (%8.1f)" if args.usecounts else " (%8.6f)"

    for n in range(1,1+args.M):
        xcent = cent[:n] ## all centroids so far
        xcen_epicount = eu.epitope_tally(xcent)
        xcov = oc.offby_count_with_obarray(seq_epicount,
                                           xcen_epicount,obarray)

        print(f"{n:3d}",end="")

        if args.B:
            xlist= np.zeros((args.offby+1,args.B))
            for b in range(args.B):
                newseqs = []
                Nseq = len(seqs)
                for k in range(Nseq):
                    ndx = random.randrange(0,Nseq)
                    newseqs.append(seqs[ndx])
                b_total,b_seq_epicount = epitope_tally_total(newseqs,args.usebadepi)

                b_xcov = oc.offby_count_with_obarray(b_seq_epicount,xcen_epicount,obarray)
                #print "xcov:",xcov
                for d in range(args.offby+1):
                    if args.usecounts:
                        xlist[d,b] = b_xcov[d]
                    else:
                        xlist[d,b] = b_xcov[d]/b_total

        for d in range(args.offby+1):
            if not args.usecounts:
                xcov[d] /= total
            print(fmt % xcov[d],end="")
            if args.B:
                print(b_fmt % np.std(xlist[d,:]),end="")

        if args.rare:
            rarest = min( [seq_epicount[epi] for epi in xcen_epicount] )
            print("%5d" % rarest,end="")
            rarecount = [epi for epi in xcen_epicount
                         if seq_epicount[epi] == rarest]
            print("%8d" % len(rarecount),end="")

        print()

if __name__ == "__main__":
    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)

