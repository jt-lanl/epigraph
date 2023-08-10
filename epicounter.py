'''make epicounter from aligned sequences'''

import re
from collections import Counter,defaultdict

import config ## defines globals; specifically EPIMER
DASH='-'
BADCHAR='x'

def mkepitopecounter(seqstrings,dash=DASH, weights=None, mer=0):
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

    EPIMER = mer or config.EPIMER

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
            if BADCHAR not in epi:
                epicount[n][epi] += weight

    return epicount
