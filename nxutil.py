'''nxutil: utilities for networkx'''

import random
import networkx as nx
from verbose import verbose as v

Begin = "==BEGIN=="
End   = "===END==="

def add_edge_safely(G,nodeleft,noderight):
    """
    add edge nodeleft->noderight without creating new nodes
    don't add an edge that connects a node to itself
    """
    if nodeleft in G and noderight in G and nodeleft != noderight:
        G.add_edge(nodeleft,noderight)
    else:
        v.vprint("Did NOT add edge: ",nodeleft,"->",noderight)
        if nodeleft not in G:
            v.vprint("   ...because",nodeleft,"not in G")
        if noderight not in G:
            v.vprint("   ...because",noderight,"not in G")

def cutbacklinks(G,thresh=0):
    """
    For a graph in which each node has a ppos (pseudo-position), cut
    all links between nodes A->B if A[ppos] - B[ppos] > thresh
    """
    print("Cutting back links with thresh=",thresh)
    backlinks=[]
    for e in G.edges():
        Appos = G.nodes[e[0]]['ppos']
        Bppos = G.nodes[e[1]]['ppos']
        delta_ppos = Appos - Bppos
        if (delta_ppos > thresh) and (e[1] != End):
            v.vprint("remove e:",e,delta_ppos,":",Appos,"-",Bppos)
            backlinks.append(e)
    print("Removed",len(backlinks),"backlinks")
    G.remove_edges_from(backlinks)
    ## replace backlink with edge to End, so as not to produce deadend
    for e in backlinks:
        G.add_edge(e[0],End)

def maxbacklinkinpath(G,path):
    maxbacklink = -1000
    nmax=0
    for n in range(len(path)-2):  ## -2 so as not to include End
        Appos = G.nodes[path[n]]['ppos']
        Bppos = G.nodes[path[n+1]]['ppos']
        backlink = Appos - Bppos
        if backlink > maxbacklink:
            maxbacklink=backlink
            nmax=n
    print("backlink: ",nmax,path[nmax],path[nmax+1],maxbacklink)
    return maxbacklink


def edgecyclelist(G,e,with_sibling=False):
    ## Given an edge that is part of a cycle,
    ## return the list of edges that make up the cycle
    edgelist=[]
    nprev = e[0]
    for n in nx.shortest_path(G,source=e[1],target=e[0]):
        if not with_sibling or edgesiblingtest(G,[nprev,n]):
            edgelist.append([nprev,n])
        nprev=n
    return edgelist

def edgepposdiff(G,e):
    ## corresponds to backlink distance
    return G.nodes[e[0]]['ppos']-G.nodes[e[1]]['ppos']

def edgepposval(G,e):
    return G.nodes[e[0]]['ppos']

def edgecountsum(G,e):
    return G.nodes[e[0]]['count']+G.nodes[e[1]]['count']
    #return G.nodes[e[1]]['count']

def edgecountmax(G,e):
    return max( [G.nodes[e[0]]['count'],G.nodes[e[1]]['count']] )

def edgesiblingtest(G,e):
    return len(G.pred[e[1]])>1 or len(G.succ[e[0]])>1

def edgespousetest(G,e):
    # edge (e0,e1) has a "spouse" if e1 has more parents than just e0
    return len(G.pred[e[1]])>1

def edgeheuristic(G,e):
    ## triple heuristic
    val = 1.0e10
    if edgesiblingtest(G,e):
    #if edgespousetest(G,e): ## was edgesiblingtest
        val =  edgecountmax(G,e) #- 1.0e-6*edgepposdiff(G,e)
        #val =  edgecountsum(G,e) #- 1.0e-6*edgepposdiff(G,e)
        #val =  G.nodes[e[0]]['count']
        #val =  1.0e-6*edgecountsum(G,e) - edgepposdiff(G,e)
    return val

def edge_new_heuristic(G,e):
    ## distinguish between spouse and sibling
    ## if spouse AND sibling, that's best, score=0
    ## if spouse but not sibling, then e[0] is isolated
    ## if sibling but not spouse, then e[1] is isolated
    ## if neither sibling nor spouse, then that's worst
    if len(G.pred[e[1]])>1:    ## if spouse
        if len(G.succ[e[0]])>1:  ## if sib
            val = 0.5*edgecountmax(G,e)
        else:
            val = edgecountmax(G,e) #G.nodes[e[0]]['count']
    else:                              ## if not spouse
        if len(G.succ[e[0]])>1:  ## if sib
            val = edgecountmax(G,e) #G.nodes[e[1]]['count']
        else:
            val = 1.0e10
    return val

def edge_ss_heuristic(G,e):
    ## distinguish between spouse and sibling
    ## if spouse AND sibling, that's best, score=0
    ## if spouse but not sibling, then e[0] is isolated
    ## if sibling but not spouse, then e[1] is isolated
    ## if neither sibling nor spouse, then that's worst
    val = edgecountsum(G,e)
    if len(G.pred[e[1]])==1:   ## if not spouse
        val += G.nodes[e[1]]['count']
    if len(G.succ[e[0]])==1:    ## if not sib
        val += G.nodes[e[0]]['count']
    return val



def min_edgeheuristic(G,e,heuristic=edgeheuristic,with_sibling=False,sign=1):
    leastval= 1.0e10
    for ee in edgecyclelist(G,e,with_sibling=with_sibling):
        val = sign*heuristic(G,ee)
        if val < leastval:
            leastval=val
            emin = ee
    return emin

def min_edgeheuristic_rnd(G,e,heuristic=edgeheuristic,with_sibling=False,sign=1):
    eelist = edgecyclelist(G,e,with_sibling=with_sibling)
    vlist = [sign*heuristic(G,ee) for ee in eelist]
    vmin = min(vlist)
    vargminlist = [i for i,v in enumerate(vlist) if v == vmin]
    vargmin = random.choice(vargminlist)
    return eelist[vargmin]

def max_edgeheuristic(G,e,heuristic=edgeheuristic,with_sibling=False):
    return min_edgeheuristic(G,e,heuristic=heuristic,
                             with_sibling=with_sibling,sign=-1)


def dagify(G,Heuristic='SNC3'):
    ## Remove cycle-inducing edges until G is a DAG
    ## Note that the "optimal" way to do this is an NP-hard problem
    ## Various heuristic approaches are employed here; make no claims of optimality

    ## Modifies G, use "Gnew = dagify(G.copy())"
    ## if you want to keep the old G

    count_removed_edges=0
    max_removed_count=0
    sum_removed_count=0
    while True:
        e = find_next_cycle_inducing_edge_alt(G)
        if not e:
            print("DAGIFY: Removed",count_removed_edges,
                  "edges with values: max=",max_removed_count,
                  "; sum=",sum_removed_count)
            return G

        if Heuristic == 'LBL':
            ## Heuristic, cut with largest ppos diff (LBL = Largest BackLink)
            e = max_edgeheuristic(G,e,heuristic=edgepposdiff)

        if Heuristic == 'LPP':
            ## Another Heuristic, cut with largest ppos (LPP = Largest PPos)
            ## This works surprisingly well, given how simple it is!
            e = max_edgeheuristic(G,e,heuristic=edgepposval)

        if Heuristic == 'LPS':
            ## Cut with largest ppos, assuming sibling exists (LPS = Largest Ppos w/ Sibling)
            e = max_edgeheuristic(G,e,heuristic=edgepposval,with_sibling=True)

        if Heuristic == 'RSNC':
            ## SNC = Random Smallest Node Count
            ## Cut edges whose node counts are smallest
            e = min_edgeheuristic_rnd(G,e,heuristic=edgecountsum)

        if Heuristic == 'RSSNC':
            ## RSSNC = Random Sibling Smallest Node Count
            ## Cut edges (that have siblings or spouses) whose node counts are smallest
            e = min_edgeheuristic_rnd(G,e,heuristic=edgecountsum,with_sibling=True)

        if Heuristic == 'FWS':
            ## Get the first edge which has a sibling (FWS = First edge With Sibling)
            for ee in edgecyclelist(G,e):
                if edgesiblingtest(G,ee):
                    e = ee
                    break

        if Heuristic == 'FWSP':
            ## Get the first edge which has a spouse (FWSP = First edge With SPouse)
            for ee in edgecyclelist(G,e):
                if edgespousetest(G,ee):
                    e = ee
                    break

        if Heuristic == 'RWS':
            ## Get a random edge which has a sibling (RWS = Random With Sibling)
            e = random.choice( edgecyclelist(G,e,with_sibling=True) )

        if Heuristic == 'RND':
            ## RND = Random
            e = random.choice(edgecyclelist(G,e))

        if Heuristic == 'SNC3':
            ## Combine three heuristics:
            ## 1. avoid isolating nodes by ensuring
            ##    multiple predecessors to e[1]
            ## 2. cut edges whose node counts are smallest
            ## 3. breaking ties according to the longest backlink
            e = min_edgeheuristic(G,e,heuristic=edgeheuristic)

        if Heuristic == 'RSC2':
            ## Combine two heuristics:
            ## 1. avoid isolating nodes by ensuring multiple predecessors to e[1] (spouse)
            ## 2. cut edges whose node counts are smallest
            ## 3. randomly choose from among the ties
            e = min_edgeheuristic_rnd(G,e,heuristic=edgeheuristic)

        if Heuristic == 'RSS':
            ## Random Spouse/Sibling
            ## Cut edges with low node counts, with extra counts
            ## given to non-spouse and non-sibling edges
            ## "Between" RSNC and RSSNC
            e = min_edgeheuristic_rnd(G,e,heuristic=edge_ss_heuristic)

        edgecountval = edgecountsum(G,e)
        v.print("Actually remove edge: ",e,end=" ")
        v.print("value:",edgecountval,end=" ")
        v.print("ppos:",G.nodes[e[0]]['ppos'],"->",G.nodes[e[1]]['ppos'])
        G.remove_edge(*e)
        count_removed_edges += 1
        sum_removed_count += edgecountval
        max_removed_count = max( [max_removed_count, edgecountval] )


def remove_cycle_inducing_edge(G):
    ## Do a topological sort, but
    ## where the sort identifies a cycle, break the cycle
    ## Return True if G is a DAG, and G will be unchanged
    ## Retrun False if an edge was removed; then it might or
    ## might not be a DAG, you'll just have to run it again to see!

    ## copied from NetworkX implementation of topological_sort
    ## and modified to remove the first cycle-inducing edge

    seen = set()
    explored = set()

    for vv in G.nodes():
        #print "vv=",vv
        if vv in explored:
            continue
        fringe = [vv]   # nodes yet to look at
        while fringe:
            #print "fringe:",fringe
            w = fringe[-1]  # depth first search
            if w in explored: # already looked down this branch
                fringe.pop()
                continue
            seen.add(w)     # mark as seen
            # Check successors for cycles and for new nodes
            #print "w:",w,
            new_nodes = []
            for n in G.successors(w): # shorthand n in G[w]:
                ## if n in explored, then we've seen n and all its successors
                if n not in explored:
                    if n in seen: #CYCLE !!
                        print("remove cycle-inducing edge: ",w,"->",n)
                        G.remove_edge(w,n)
                        ## If n has no predecessors, give it Begin
                        if not G.predecessors(n):
                            G.add_edge(Begin,n,weight=0)
                        return False
                        #raise nx.NetworkXUnfeasible("Graph contains a cycle.")
                    new_nodes.append(n)
            #print " -> new nodes:",new_nodes
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored.add(w)
                fringe.pop()    # done considering this node
    return True

def find_cycle_inducing_edge(G):
    ## Do a topological sort, but
    ## where the sort identifies a cycle, break the cycle
    ## Return True if G is a DAG, and G will be unchanged
    ## Return the edge (w,n) that induces a cycle

    ## copied from NetworkX implementation of topological_sort
    ## and modified to return the first cycle-inducing edge

    seen = set()
    explored = set()

    for vv in G.nodes():
        if vv in explored:
            continue
        fringe = [vv]   # nodes yet to look at
        while fringe:
            w = fringe[-1]  # depth first search
            if w in explored: # already looked down this branch
                fringe.pop()
                continue
            seen.add(w)     # mark as seen
            # Check successors for cycles and for new nodes
            new_nodes = []
            for n in G.successors(w): # shorthand n in G[w]:
                ## if n in explored, then we've seen n and all its successors
                if n not in explored:
                    if n in seen: #CYCLE !!
                        return (w,n)
                        #raise nx.NetworkXUnfeasible("Graph contains a cycle.")
                    new_nodes.append(n)
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored.add(w)
                fringe.pop()    # done considering this node
    return None

def find_cycle_inducing_edge_rnd(G):
    ## Do a topological sort, but
    ## where the sort identifies a cycle, break the cycle
    ## Return True if G is a DAG, and G will be unchanged
    ## Return the edge (w,n) that induces a cycle

    ## copied from NetworkX implementation of topological_sort
    ## and modified to return the first cycle-inducing edge

    seen = set()
    explored = set()

    nodelist = G.nodes()
    random.shuffle(nodelist)

    for vv in nodelist: #G.nodes():
        if vv in explored:
            continue
        fringe = [vv]   # nodes yet to look at
        while fringe:
            w = fringe[-1]  # depth first search
            if w in explored: # already looked down this branch
                fringe.pop()
                continue
            seen.add(w)     # mark as seen
            # Check successors for cycles and for new nodes
            new_nodes = []
            for n in G.successors(w): # shorthand n in G[w]:
                ## if n in explored, then we've seen n and all its successors
                if n not in explored:
                    if n in seen: #CYCLE !!
                        return (w,n)
                        #raise nx.NetworkXUnfeasible("Graph contains a cycle.")
                    new_nodes.append(n)
            if new_nodes:   # Add new_nodes to fringe
                fringe.extend(new_nodes)
            else:           # No new nodes so w is fully explored
                explored.add(w)
                fringe.pop()    # done considering this node
    return None



global_clist=[] ## ack! global
global_ccc=0    ## count calls to scc
def find_next_cycle_inducing_edge_alt(G):
    ## Use Tarjan's algorithm
    global global_ccc
    if not global_clist:
        global_ccc += 1
        for comp in nx.strongly_connected_components(G):
            if len(comp)>1:
                global_clist.append(comp)
    if global_clist:
        comp = global_clist.pop()
        #global_clist[:] = global_clist[1:]
        #print global_ccc,"edge:",comp[0:2],"of length",len(comp),
        #"; still ",len(global_clist),"remaining nontrivial components"
        ## Question: is the component itself a cycle?? is it in cycle order?
        (a,b) = comp[0:2]
        p = nx.shortest_path(G,source=a,target=b)
        #print "edge:",p[0:2],"of length",len(comp),len(p),
        #len(nx.shortest_path(G,source=p[1],target=p[0]))
        return tuple(p[0:2])
    else:
        return None


def find_cycle_inducing_edge_alt(G):
    ## Use Tarjan's algorithm (fresh each time, but only take the first one)
    for comp in nx.strongly_connected_components(G):
        if len(comp)>1:
            ## Choose nodes two at random ... amazing how badly this works
            ## (score ok, but runtime is impossible!)
            #(a,b) = random.sample(comp,2)
            (a,b) = comp[0:2]
            p = nx.shortest_path(G,source=a,target=b)
            print("edge:",p[0:2],"of length",len(comp),len(p),
                  len(nx.shortest_path(G,source=p[1],target=p[0])))
            return tuple(p[0:2])
    return None ## does this ever happen?



## From Guido's personal blog!
## http://neopythonic.blogspot.com/2009/01/detecting-cycles-in-directed-graph.html
## A little simpler than using the topological sort (though same basic idea)
## And it returns the whole cycle of nodes, not just the cycle-inducing edge
def find_cycle_guido(G):
    todo = set(G.nodes())
    while todo:
        node = todo.pop()
        stack = [node]
        while stack:
            top = stack[-1]
            ## node within node??
            for node in G.successors(top):
                if node in stack:
                    return stack[stack.index(node):]
                if node in todo:
                    stack.append(node)
                    todo.remove(node)
                    break
            else:
                node = stack.pop()
    return None


## Another web-obtained routine (probably equivalent to nx's version)
def strongly_connected_components(graph):
    """
    Tarjan's Algorithm (named for its discoverer, Robert Tarjan)
    is a graph theory algorithm for finding the strongly connected
    components of a graph.

    Based on: http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    """

    index_counter = [0]
    stack = []
    lowlinks = {}
    index = {}
    result = []

    def strongconnect(node):
        # set the depth index for this node to the smallest unused index
        index[node] = index_counter[0]
        lowlinks[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)

        # Consider successors of `node`
        try:
            successors = graph[node]
        except:
            successors = []
        for successor in successors:
            if successor not in lowlinks:
                # Successor has not yet been visited; recurse on it
                strongconnect(successor)
                lowlinks[node] = min(lowlinks[node],lowlinks[successor])
            elif successor in stack:
                # the successor is in the stack and hence
                # in the current strongly connected component (SCC)
                lowlinks[node] = min(lowlinks[node],index[successor])

        # If `node` is a root node, pop the stack and generate an SCC
        if lowlinks[node] == index[node]:
            connected_component = []

            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            component = tuple(connected_component)
            # storing the result
            result.append(component)

    for node in graph:
        if node not in lowlinks:
            strongconnect(node)

    return result


def clear_node_max_values(G):
    for n in G.nodes():
        G.nodes[n]['value']=0
        G.nodes[n].pop('pre',None)

def maximum_node_bandwidth(G,topnodes=False):
    """
    Find the largest value N such that a path exists through the graph
    with G[node]['count'] >= N for all nodes in the path
    """

    if not topnodes:
        topnodes = nx.topological_sort(G)

    bandwidth=dict()
    bandwidth[Begin]=999999 ## Infinity

    for node in topnodes:
        if node == Begin:
            continue

        maxpredval=max(bandwidth[p] for p in G.predecessors(node))

        if node == End:
            bandwidth[End] = maxpredval
        else:
            bandwidth[node] = min([maxpredval,
                                   G.nodes[node]['count']
            ])

    return bandwidth[End]


def ok_populate_node_max_values(G,topnodes=False):

    clear_node_max_values(G) ## for safety's sake

    if topnodes is False:
        v.print("pop-max sorting...")
        topnodes = nx.topological_sort(G)

    G.nodes[Begin]['value']=0
    G.nodes[Begin]['pre']=[]

    for currentnode in topnodes:
        if currentnode == Begin:
            continue

        ## find best predecessors to the current node

### This is a more concise (but, oddly, more expensive) way to do the same thing
#        pnodes = list(G.predecessors(currentnode))
#        cvals = [G.nodes[pnode]['value'] for pnode in pnodes]
#        maxval = max(cvals)
#        bestpred = [pnodes[n] for n in range(len(cvals)) if cvals[n]==maxval]
### How about?
#        maxval = max(G.nodes[p]['value'] for p in G.precessors(currentnode))
#        bestpred = [p for p in G.pred... if G.nodes[p][value] == maxval]

        maxval=-1
        bestpred = []
        for pnode in G.predecessors(currentnode):

            cval = G.nodes[pnode]['value']
            if cval > maxval:
                maxval = cval
                bestpred = []
            if cval >= maxval:
                bestpred.append(pnode)

        G.nodes[currentnode]['value']=maxval + G.nodes[currentnode]['count']
        G.nodes[currentnode]['pre']=bestpred ## list of best predecessors
        #print "[",currentnode,"]:",maxval,", via",bestpred

def populate_node_max_values(G,topnodes=False):

    clear_node_max_values(G) ## for safety's sake

    if topnodes is False:
        v.print("pop-max sorting...")
        topnodes = nx.topological_sort(G)

    G.nodes[Begin]['value']=0
    G.nodes[Begin]['pre']=[]

    for currentnode in topnodes:
        if currentnode == Begin:
            continue

        ## find best predecessors to the current node
        pred = list(G.predecessors(currentnode))
        maxval=max(G.nodes[p]['value'] for p in pred)
        bestpred = [p for p in pred if G.nodes[p]['value'] == maxval]

        G.nodes[currentnode]['value']=maxval + G.nodes[currentnode]['count']
        G.nodes[currentnode]['pre']=bestpred ## list of best predecessors


def optimal_path_graph_single(G):
    ## like optimal_path_graph, but with a single path
    path = []
    node = End
    while True:
        path.append(node)
        if node == Begin:
            break
        nodelist = G.nodes[node]['pre']
        node = nodelist[0]

    H = nx.DiGraph()
    prevnode=None
    for node in path:
        H.add_node(node)
        if prevnode:
            H.add_edge(node,prevnode)
        prevnode=node
    return H


def optimal_path_graph(G):
    """
    Returns a graph H that is a subgraph of G with the property that
    /all/ paths through H are /optimal/ paths through G
    """
    H = nx.DiGraph()
    CurrentNodes=set() # use set, not list, avoid duplicated nodes
    CurrentNodes.add(End)
    while CurrentNodes:
        node = CurrentNodes.pop()
        H.add_node(node)
        ## Is this the same as G.pred[node] ? No, only want best pred's !
        for prenode in G.nodes[node]['pre']:
            H.add_edge(prenode,node)
            CurrentNodes.add(prenode)
    return H

def optimal_path_graph2(G):
    """
    Returns a graph H which is a subset of G with the property that
    /all/ paths through H are /optimal/ paths through G
    """
    ## This version might be faster than the simpler optimal_path_graph; but
    ## now that we are using set instead of list, the simple is fast enough
    ## Do this in two steps: first identify the nodes in a set, then make graph
    Hnodeset=set()
    CurrentNodes=[End]
    while CurrentNodes:
        node = CurrentNodes.pop()
        Hnodeset.add(node)
        for prenode in G.nodes[node]['pre']:
            if prenode not in Hnodeset:
                CurrentNodes.append(prenode)

    H = nx.DiGraph()
    for node in Hnodeset:
        H.add_node(node)
        for prenode in G.nodes[node]['pre']:
            H.add_edge(prenode,node)

    return H

def find_an_optimal_path(G,topnodes=False):
    """
    Return a path (list of nodes [Begin,...,End])
    for which the sum of node counts is maximal.
    There may be more than one; but just return one.
    """

    populate_node_max_values(G,topnodes)#False)
    H = optimal_path_graph(G)
    path = one_path(H)
    return path

def find_a_greedy_path(G,topnodes=False):
    """
    Return a path (list of nodes [Begin,...,End])
    that greedily takes the highest scoring successor
    """

    populate_node_max_values(G,topnodes)
    node = Begin
    path = [node]
    while node != End:
        node = max(G.sucessors(node),
                   key=lambda n: G.nodes[n]['count'])
        path.append(node)

    return path


def one_path(G):
    path=[]
    n = Begin
    while n != End:
        path.append(n)
        n = next(G.successors(n))  ## or maybe G.succ[n][0]
    path.append(End)
    return path


def clear_ppos(G):
    for n in G.nodes():
        G.nodes[n].pop('ppos',None)

def set_ppos(G):
    """
    ppos is pseudo-position; defined by
    the number of nodes between it and the Begin node
    """
    nodelist=[Begin]
    ppos=0
    for n in nodelist:
        G.nodes[n]['ppos']=ppos
    while nodelist:
        ppos +=1
        newlist=[]
        for n_cur in nodelist:
            for n in G.successors(n_cur):
                if 'ppos' not in G.nodes[n]:
                    G.nodes[n]['ppos'] = ppos ## get first value
                    newlist.append(n)          ## only append if new
        nodelist=newlist

def ppos_isolated(G):
    ## even after removing deadends and deadbegins,
    ## there may still be isolated cycles; nodes in
    ## these isolated cycles will lack a 'ppos' key
    return [n for n in G.nodes()
            if 'ppos' not in G.nodes[n].keys()]

def remove_ppos_isolated(G):
    isolated = ppos_isolated(G)
    v.print("Removing",len(isolated),"isolated nodes")
    for n in isolated:
        G.remove_node(n)
    return G

def list_deadends(G,among_nodes=None):
    if among_nodes is None:
        among_nodes = G.nodes()
    deadlist = [epi for epi in among_nodes
                if epi!=End
                and not G.succ[epi]] #list(G.successors(epi))]
    return deadlist

def list_deadbegins(G,among_nodes=None):
    if among_nodes is None:
        among_nodes = G.nodes()
    deadlist = [epi for epi in among_nodes
                if epi != Begin
                and not G.pred[epi]] #list(G.predecessors(epi))]
    return deadlist

def remove_deadends(G):
    deadlist = list_deadends(G)
    while deadlist:
        pre_dead = set()
        for n in deadlist:
            pre_dead.update(G.predecessors(n))
            G.remove_node(n)
        deadlist = list_deadends(G,among_nodes=pre_dead)
    return G

def remove_deadbegins(G):
    deadlist = list_deadbegins(G)
    while deadlist:
        post_dead = set()
        for n in deadlist:
            post_dead.update(G.successors(n))
            G.remove_node(n)
        deadlist = list_deadbegins(G,among_nodes=post_dead)
    return G


def random_path(G):
    return random_path_to_end(G,Begin)

def random_path_to_end(G,startnode):
    p=[]
    node=startnode
    p.append(node)
    while node != End:
        successor_nodes = list( G.successors(node) )
        if not successor_nodes:
            raise RuntimeError(f"Dead End -- No successors to node {node}"
                               f"with count={G.nodes[node]['count']}\n"
                               f"p={p}")
        node = random.choice( successor_nodes )
        p.append(node)
    return p

def random_path_from_begin(G,endnode):
    p=[]
    node=endnode
    p.insert(0,node)
    while node != Begin:
        predecessor_nodes = list(G.predecessors(node))
        if not predecessor_nodes:
            raise RuntimeError(f"Dead End -- No predecessors to node {node}"
                               f"with count={G.nodes[node]['count']}\n"
                               f"p={p}")
        node = random.choice( predecessor_nodes )
        p.insert(0,node) # prepend
    return p
