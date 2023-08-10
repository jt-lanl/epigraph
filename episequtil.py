''' sequence utilities '''

import re
from collections import Counter

import verbose as v
import sequtil


import config
import epiuutil as eu

def init_epimer_length(epimer_length):
    config.EPIMER = epimer_length

def too_many_repeats_in_seq(seq,tmr=0):
    ''' return True if too many epitopes appear more than once in the sequence '''
    if tmr==0:
        return False
    episeq = eu.episeq_from_string(seq,rmdash=True)
    repeats_in_str = sum( (episeq.count(epi)>1) for epi in set(episeq) )
    ## but not max(episeq.count(epi) for epi in set(episeq))
    ## or: any(episeq.count(epi) >= tmr for epi in set(episeq))
    ## which would be if any individual epitopes appear too many times
    return bool(repeats_in_str >= tmr)


def seqgen_filter(seqs,reverse=False,tmr=0,R=False,rmx=True):
    ''' filter incoming sequences according to various criteria '''

    if R and len(R)!=2:
        R=False
    if R and len(R)==2:
        lo,hi = R
        assert lo < hi

    for s in seqs:
        if rmx or tmr:
            e = eu.episeq_from_string(s.seq,rmdash=True)
            if rmx and len(e)==0:
                v.print("No epitopes, so ignoring seq:",s.name)
                continue
            if tmr:
                nrepeats = sum( e.count(epi)>1 for epi in set(e) )
                if nrepeats > tmr:
                    v.print("Too many repeats",nrepeats,"in seq:",s.name)
                    v.print([(epi,e.count(epi)) for epi in set(e) if e.count(epi)>1])
                    continue
        if R:
            if lo > len(s.seq):
                raise RuntimeError(f"ERROR! lo too high: {lo} > {len(s.seq)}")
            s.seq = s.seq[lo:hi]

        if reverse:
            s.seq = s.seq[::-1]
        yield s

def seq_filter(seqs,reverse=False,tmr=0,R=0,rmx=True):
    '''
    reverse=False: reverse all the strings
    tmr=0: filter out strings with too many repeats
    R=0: if R=(lo,hi) then truncate strings to this range
    rmx: filter out strings that have no epitopes
    '''
    if reverse:
        for s in seqs:
            s.seq = s.seq[::-1]
    if tmr:
        seqs = eu.seq_filter_repeats(seqs,tmr)

    if R and len(R)==2:
        v.print("Truncate to range:",R)
        seqs = truncate_strings(seqs,R)

    if rmx:
        badseqs=[]
        for s in seqs:
            e = eu.episeq_from_string(s.seq,rmx=True)
            if len(e)==0:
                v.print("No epitopes, so ignoring seq:",s.name)
                badseqs.append(s)
        if badseqs:
            print("Ignoring",len(badseqs),"bad sequences")
            seqs = [s for s in seqs if s not in badseqs]

    return seqs

def seq_report_stats(seqs):
    '''
    print out various statistics on the sequences
    '''
    print("Read",len(seqs),"sequences")
    seqlen = [len(s.seq) for s in seqs]
    print("Sequences have between",min(seqlen),"and",max(seqlen),"characters")

    ## Histogram of sequence lengths
    seqlend=Counter()
    for s in seqs:
        lenstr = len(s.seq)
        seqlend[lenstr] += 1
    for lenstr in sorted(seqlend.keys()):
        print("len = %d, number of sequences: %d"%(lenstr,seqlend[lenstr]))


## Write cocktail of paths to file args.out
def cocktail_to_file(cocktail,filename,seqnames=None,reverse=False):
    if not filename:
        return
    v.print("Writing cocktail of",len(cocktail),"sequences to file:",filename)

    seqs = []
    for m,e in enumerate(cocktail):
        ## A cocktail is composed of episeq's (sequence of epitopes)
        s = eu.episeq_to_string(eu.path_chomp(e))
        if reverse:
            s = s[::-1]
        if seqnames and len(seqnames) > m:
            name = seqnames[m]
        else:
            name = "EGN-%d"%(m+1)
        seqs.append( sequtil.SequenceSample( name, s ) )
    sequtil.write_seqfile(filename,seqs)

def check_duplicated_epi(s):
    dups=set()
    pcount = Counter()
    p = eu.episeq_from_string(s,rmdash=True)
    for epi in p:
        pcount[epi] += 1
        if pcount[epi]>1:
            print("Duplicated epitope: ",epi)
            dups.add(epi)
    return dups



def remove_bad_epitopes(epicount):
    for epi in list(epicount.keys()):
        if 'x' in epi:
            epicount.pop(epi)
        if 'X' in epi:
            raise RuntimeError("Where did that X come from?? " + epi)

def truncate_strings(seqs,lohi):
    lo,hi = lohi
    if lo>hi:
        raise RuntimeError("ERROR: lo>hi")
    for seq in seqs:
        if lo > len(seq.seq):
            raise RuntimeError("ERROR! lo too high: %d > %d" % (lo,len(seq.seq)))
        seq.seq = seq.seq[lo:hi]
    return seqs



def make_posn_adjust_map(seqstr):
    """ returns an array posn_adjust[] such that
    long_offset = posn_adjust[short_offset]
    where long_offset is position in the original string (seqstr)
    and short_offset is position in dash-removed string
    """
    return [ n_long for n_long,seqchar in enumerate(seqstr)
             if seqchar != "-" ]

def make_epi_posn_dict(long_seqstr):
    EPIMER = config.EPIMER
    re_dash=re.compile('-')
    short_seqstr = re_dash.sub('',long_seqstr)
    posn_adjust = make_posn_adjust_map(long_seqstr)

    ## dictionary of positions of epitopes in the long string
    ## (what if same epitope appears in two places??)
    posn_in_long_string_of_epi={}

    for off_short in range(len(short_seqstr)-EPIMER+1):
        off_long = posn_adjust[off_short]
        end_long = posn_adjust[off_short+EPIMER-1]
        epi = long_seqstr[off_long:end_long+1]
        if "x" in epi:
            continue
        posn_in_long_string_of_epi[epi] = off_long

    return posn_in_long_string_of_epi



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",help="input tbl filename")
    args = parser.parse_args()

    seqs = read_tbl(args.filename)

    print("Read",len(seqs),"sequences")
    print("Name of first sequence:",seqs[0].name)
    print("Name of last sequence: ",seqs[-1].name)
    print("Sequence length:",len(seqs[0].seq))
    print("Sequence length:",len(seqs[1].seq))
    print("Sequence length:",len(seqs[2].seq))
    print("Sequence length:",len(seqs[-1].seq))
    print("Last sequence:",seqs[-1].seq)


    posn_adj = make_posn_adjust_map(seqs[0].seq)
    print(posn_adj)

    import pylab
    pylab.plot(posn_adj)
    pylab.show()
