'''Tabulate positions for each kmer that appears in an alignment'''

from collections import defaultdict
from epicounter import mkepitopecounter

class MerPositionTable:
    '''make table of position (actually, list of positions) for each K-mer'''

    def __init__(self,K=9):
        self.K=K
        self.table = defaultdict(list)

    def update_with_epicount(self,epicount):
        '''update the table based on an epicount structure'''
        ## epicount[pos][epi] is count of seqs where epi appears at pos
        for pos in epicount:
            for epi in epicount[pos]:
                self.table[epi].append(pos)
        return self

    def update_with_seqstrings(self,seqstrings):
        '''update the table based on a list of sequence strings'''
        ## nb: seqstrings = [s.seq for s in seqlist]
        epicount = mkepitopecounter(seqstrings,mer=self.K)
        return self.update_with_epicount(epicount)

    def update_with_seqlist(self,seqlist):
        '''update the table based on a list of sequences'''
        seqstrings = [s.seq for s in seqlist]
        return self.update_with_seqstrings(seqstrings)

    def get_positions_for_epi(self,epi):
        '''return list of positions at which epi occurs'''
        return self.table[epi]

    def get_positiosn_for_epilist(self,epilist):
        '''return list of positions covered by all epi's'''
        positions = set()
        for epi in epilist:
            positions.update(self.get_positions_for_epi(epi))
        return sorted(positions)
    
