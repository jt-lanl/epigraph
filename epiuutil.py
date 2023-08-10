''' epiuutil: EPI-graph, U-naligned, UTIL-ities'''

import re
from collections import Counter

import networkx as nx
import config
from verbose import verbose as v
import decycle

import nxutil
Begin = nxutil.Begin
End   = nxutil.End

def episeq_from_string(seqstr,rmx=True,rmdash=True):
    '''
    input a seq string and get a list of epitopes of length config.EPIMER
    rmdash: remove '-' from seqstring before creating list of epitopes
    rmx: remove epitopes that contain 'x'
    '''
    EPIMER=config.EPIMER
    if rmdash:
        re_dash=re.compile('-')
        seqstr = re_dash.sub('',seqstr)
    episeq = [seqstr[n:n+EPIMER] for n in range(len(seqstr)-EPIMER+1)]
    if rmx and 'x' in seqstr:
        episeq = [e for e in episeq if 'x' not in e]
    return episeq

def episeq_to_string(episeq):
    '''
    Convert a consistent list [ABC, -BC, BCD, CDE ...] to string "A-BC..."
    '''
    seqstr = "".join(epi[0] for epi in episeq) ## first char of each epitope
    seqstr += episeq[-1][1:] ## all but first char [1:] of last epitope [-1]
    return seqstr

def epitope_tally(seqs,rmx=True,rmdash=True):
    ''' Input is a list of SequenceSamples
    Returns a Counter object, epicount,
    for which epicount[epitope] is the number of
    sequences in which the given epitope appears
    '''
    epicount = Counter()
    for s in seqs:
        episeq = episeq_from_string(s.seq,rmx=rmx,rmdash=rmdash)
        epicount.update(set(episeq))
    return epicount

def seq_filter_repeats(seqs,tmr):
    ''' Return a list of SequenceSample's which avoids those
    sequences with too many repeats (more than tmr)
    '''

    if tmr == 0:
        return seqs

    sslist = []
    for s in seqs:
        episeq = episeq_from_string(s.seq,rmdash=True)
        repeats = sum(1
                      for epi in set(episeq)
                      if episeq.count(epi) > 1)
        if repeats < tmr:
            sslist.append(s)
    if len(seqs) > len(sslist):
        v.vprint(f"Read {len(seqs)} sequences, "
               f"Removed {len(seqs)-len(sslist)} sequences, "
               f"Kept {len(sslist)} sequences.")
    return sslist

def graph_cocktail_coverage(G,c):
    epi_covered = set( epi for path in c for epi in path )
    epi_covered &= set(G.nodes) ## in case there are rogue nodes in c
    return sum( G.nodes[epi]['count'] for epi in epi_covered )

def cocktail_coverage(c,epicount):
    epi_covered=set()
    for p in c:
        for epi in p:
            epi_covered.add(epi)
    return sum( [epicount.get(epi,0) for epi in epi_covered] )

def path_chomp(path):
    ''' Given a path, remove Begin and End nodes if they are present
    '''
    if path[0] == nxutil.Begin:
        path = path[1:]
    if path[-1] == nxutil.End:
        path = path[:-1]
    return path

def remove_lowcount_nodes(epicount,lowcount):
    '''
    Exclude rare epitopes, ie
    Remove nodes whose count is <= lowcount
    '''
    if lowcount==0:
        return epicount
    epinew=Counter() ## more efficient to pop bad nodes?
    for epi in epicount:
        if epicount[epi] > lowcount:
            epinew[epi] = epicount[epi]
    if len(epicount) > len(epinew):
        print("Removed",len(epicount)-len(epinew),"nodes with count <=",lowcount)
    return epinew

def remove_weak_nodes(G):
    ''' Remove weak nodes; a node is weak if
    1. It has 'epicount' value of 1
    2. It has a single   successor, and that   successor has epicount value of 1
    3. It has a single predecessor, and that predecessor has epicount value of 1
    '''
    listweaknodes=[]
    for n in G:
        if G.nodes[n]['count']>1:
            continue
        if len(G.successors(n))!=1 or len(G.predecessors(n))!=1:
            continue
        k=G.successors(n)[0]
        if G.nodes[k]['count']>1:
            continue
        k=G.predecessors(n)[0]
        if G.nodes[k]['count']>1:
            continue
        if n in (Begin,End):
            continue
        ## if we make it this far,
        ## then n is a weak node
        listweaknodes.append(n)


    if listweaknodes:
        print("Removing",len(listweaknodes),"weak nodes")
        for n in listweaknodes:
            G.remove_node(n)
    return listweaknodes

def restore_weak_nodes(G,listweaknodes,epicount):
    if not listweaknodes:
        return
    print("Putting back the weak nodes")
    for n in listweaknodes:
        G.add_node(n,count=epicount[n])

def remove_lowcount_edges(G,LOWCOUNTEDGE,epicount):
    '''
    Remove edges that connect nodes, both of whose counts are low
    '''
    listweakedges=[]
    for e in G.edges_iter():
        if epicount[e[0]] <= LOWCOUNTEDGE and \
           epicount[e[1]] <= LOWCOUNTEDGE:
            listweakedges.append(e)
    print("Remove",len(listweakedges),"weak edges")
    for e in listweakedges:
        G.remove_edge(*e)
    return listweakedges

def restore_lowcount_edges(G,listweakedges):
    if not listweakedges:
        return
    print("Putting back the weak edges")
    for e in listweakedges:
        nxutil.add_edge_safely(G,*e)
        ## Hmm, maybe we should make sure nodes exist? (and add them back too)

def connect_nodes(G,epi_right_hash=None):
    ## Our next task is to take all the epitopes, which are nodes,
    ## In fact, to take all the pairs of nodes, and connect the consistent ones
    ## Two epitopes epi_left and epi_right are consistent (and connected) if
    ## epi_left[1:] == epi_right[:-1]

    episet = set() # = set(G.nodes()) - set([Begin,End])
    for n in G:
        if n not in (Begin,End):
            episet.add(n)

    if 0:
        ## Ok, this is kind of expensive, but...
        ## It explains what we are trying to do
        for epi_left in episet:
            for epi_right in episet:
                if epi_left[1:] == epi_right[:-1]:
                    nxutil.add_edge_safely(G,epi_left,epi_right)
                    #print epi_left,"->",epi_right

    ## This may be a little confusing but it's a lot faster
    ## It uses hashing (specifically a dictionary of lists)
    ## to find the epi_right's associated with each epi_left
    if not epi_right_hash:
        epi_right_hash=dict()
        for epi in episet:
            ## all possible keys are initialized with empty list
            epi_right_hash[epi[1:]]=[]
            epi_right_hash[epi[:-1]]=[]
        for epi in episet:
            epi_right_hash[epi[:-1]].append(epi)

    ## Add edges to the graph if epi_left[1:] == epi_right[:-1]
    for epi_left in episet:
        for epi_right in epi_right_hash[epi_left[1:]]:
            nxutil.add_edge_safely(G,epi_left,epi_right)

    return epi_right_hash

def make_graph(seqs,dagify='rss',
               weights=None,
               nodelowcount=0,edgelowcount=0,rmweaknodes=False):

    ## Begin compiling epitopes
    epicount = Counter()   ## maybe defaultdict(float) instead for weights??
    firstepi = set()
    lastepi = set()

    seq_count=0
    for s in seqs:
        episeq = episeq_from_string(s.seq,rmx=True)
        if len(episeq)==0:
            v.vprint("bad sequence:",s.name)
            continue
        firstepi.add(episeq[0])
        lastepi.add(episeq[-1])
        seq_count += 1
        episeq = set(episeq)
        if weights:
            for epi in episeq:
                epicount[epi] += weights.get(s.name,1.0)
        else:
            epicount.update(episeq)

    v.vprint(f"Using {seq_count} sequences (after filtering) "
           f"to make {len(epicount)} epitope nodes")

    epi_total = sum( epicount.values() )

    ## Only keep nodes if count > NODE_LOW_COUNT
    epicount = remove_lowcount_nodes(epicount,nodelowcount)

    ## Make the (directed, possibly acyclic) Graph

    G = graph_from_epicount(epicount,firstepi,lastepi)

    G.graph['epiTotal']=epi_total

    listweaknodes=[]
    if rmweaknodes:
        listweaknodes = remove_weak_nodes(G)

    listweakedges=[]
    if edgelowcount:
        listweakedges = remove_lowcount_edges(G,edgelowcount,epicount)

    ## If no path through G, then put back those weak edges
    if not nx.has_path(G,Begin,End):
        restore_lowcount_edges(G,listweakedges)

    ## If still no path through G, then put back those weak nodes
    if not nx.has_path(G,Begin,End):
        restore_weak_nodes(G,listweaknodes,epicount)
        connect_nodes(G)
        for epi in firstepi:
            nxutil.add_edge_safely(G,Begin,epi)
        for epi in lastepi:
            nxutil.add_edge_safely(G,epi,End)

    ## If STILL no path, then something is broken
    if not nx.has_path(G,Begin,End):
        print("nodelowcount:",nodelowcount)
        print("edgelowcount:",edgelowcount)
        print("rmweaknodes:",rmweaknodes)
        with open("BadGraph.txt","w") as fout:
            print("Begin:",Begin,"->",",".join(G.successors(Begin)),file=fout)
            for node in G.nodes:
                s = G.successors(node)
                print(node,"->",",".join(s),file=fout)
        raise RuntimeError("Graph is still broken! (see BadGraph.txt)")

    G = clean_and_decycle(G,dagify=dagify)

    v.vprint(f"Final graph has {len(G)-2} nodes") ## -2 subtracts off Begin,End

    return G

def graph_from_epicount(epicount,firstepi=None,lastepi=None):

    G = nx.DiGraph()
    for epi in epicount:
        G.add_node(epi,count=epicount[epi])
    connect_nodes(G)

    ## Create and Connect to the Begin and End nodes
    G.add_node(Begin,count=0)
    if not firstepi:
        firstepi = nxutil.list_deadbegins(G)
    for epi in firstepi:
        if epi in G:
            nxutil.add_edge_safely(G,Begin,epi)

    G.add_node(End,count=0)
    if not lastepi:
        lastepi = nxutil.list_deadends(G)
    for epi in lastepi:
        if epi in G:
            nxutil.add_edge_safely(G,epi,End)

    return G


def clean_and_decycle(G,dagify):

    if not nx.has_path(G,Begin,End):
        raise RuntimeError("Graph is broken!")

    v.vvprint(f"BEFORE: G has nodes: {len(G.nodes())}")

    ## Remove those deadheads
    G = nxutil.remove_deadends(G)
    G = nxutil.remove_deadbegins(G)

    v.vvprint(f"BEFORE: G has nodes: {len(G.nodes())}")

    # Set 'ppos' = distance to Begin node (akin to position in sequence)
    nxutil.clear_ppos(G)
    nxutil.set_ppos(G)

    ## Remove isolated cycles
    G = nxutil.remove_ppos_isolated(G)

    v.vvprint(f"BEFORE: G has nodes: {len(G.nodes())}")

    ## If G has cycles, remove cycle-inducing edges, until it is a DAG
    #OLD: G = nxutil.dagify(G,args.dag)
    #NEW:
    G = decycle.decycle(G,dagify)

    v.vvprint(f"AFTER: G has nodes: {len(G.nodes())}")

    G = nxutil.remove_deadends(G) ## in case dagify introduced some dead ends
    G = nxutil.remove_deadbegins(G)

    v.vvprint(f"AFTER: G has nodes: {len(G.nodes())}")

    nxutil.clear_ppos(G)
    nxutil.set_ppos(G)
    G = nxutil.remove_ppos_isolated(G) ## shouldn't be any, but just in case

    v.vvprint(f"AFTER: G has nodes: {len(G.nodes())}")

    return G

def optimal_complementary_path(G,G_topnodes,cocktail):
    ## Set nodes in current paths to zero
    ## use restore=dict() to keep track, so those nodes can be restored

    nodes_in_current_paths = set(node
                                 for path in cocktail
                                 for node in path)

    restore=dict()
    for node in nodes_in_current_paths:
        if node in G:
            restore[node]=G.nodes[node]['count']
            G.nodes[node]['count']=0

    ## Find optimal path through this graph with zeros
    newpath = nxutil.find_an_optimal_path(G,G_topnodes)

    ## Restore nodes that were set to zero
    for node in restore:
        G.nodes[node]['count']=restore[node]

    return newpath
