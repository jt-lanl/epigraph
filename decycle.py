'''decycle: algorithms for removing cycles in a directed graph'''

import random
import time
import networkx as nx
from verbose import verbose as v
from newcycle import my_bidirectional_shortest_path

## Various edgecost functions correspond to the "cost"
## of removing the given edge

def nodecost(G,eab):
    (ea,eb)=eab
    return G.nodes[ea]['count'],G.nodes[eb]['count']

def iso_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    val = 0
    if len(G.pred[ea])==1:  ## if not spouse
        val += ca
    if len(G.succ[eb])==1:    ## if not sib
        val += cb
    return val

def supiso_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    return max((ca,cb)) + iso_edgecost(G,(ea,eb))

def rss_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    val = ca + cb
    if len(G.pred[eb])==1:  ## if not spouse
        val += cb
    if len(G.succ[ea])==1:    ## if not sib
        val += ca
    return val

def rss2_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    val = ca + cb
    if len(G.pred(eb))==1:  ## if not spouse
        val += 2*cb
    if len(G.succ[ea])==1:    ## if not sib
        val += 2*ca
    return val

def sum_sps_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    val = ca + cb
    if len(G.pred[eb])==1:  ## if not spouse
        val += cb
    return val

def sum_sib_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    val = ca + cb
    if len(G.succ[ea])==1:    ## if not sib
        val += ca
    return val



def sum_edgecost(G,eab):
    (ea,eb)=eab
    ca,cb = nodecost(G,(ea,eb))
    return ca+cb

def sup_edgecost(G,e):
    return max([G.nodes[e[0]]['count'],G.nodes[e[1]]['count']])

def rnd_edgecost(G,e):
    return 0

def pos_edgecost(G,e):
    if 'ppos' not in G.nodes[e[0]].keys():
        print("Node",e[0],"keys=",G.nodes[e[0]].keys(),"count=",G.nodes[e[0]]['count'])
        return 0
    return -1*G.nodes[e[0]]['ppos']

def dpos_edgecost(G,e):
    if 'ppos' not in G.nodes[e[0]].keys():
        return 0
    if 'ppos' not in G.nodes[e[1]].keys():
        return 0
    return G.nodes[e[1]]['ppos']-G.nodes[e[0]]['ppos']

edgecost_fcn_from_string = {
    'rss' : rss_edgecost,
    'rs2' : rss2_edgecost,
    'sum' : sum_edgecost,
    'sup' : sup_edgecost,
    'iso' : iso_edgecost,
    'supiso' : supiso_edgecost,
    'sso' : supiso_edgecost,
    'sib' : sum_sib_edgecost,
    'sps' : sum_sps_edgecost,
    'rnd' : rnd_edgecost,
    'pos' : pos_edgecost,
    'dpos' : dpos_edgecost,
}

def set_edgecost_all(G,edgecostfcn=rss_edgecost):
    ## use edge attributes to compute each cost once for each edge
    for e in G.edges:
        G.adj[e[0]][e[1]]['cost'] = edgecostfcn(G,e)

def mincost_edge_in_cycle(G, edge_list, edgecostfcn=rss_edgecost, sign=1):
#    cost_list = [sign*edgecostfcn(G,e) for e in edge_list]
    cost_list = [sign*G.adj[e[0]][e[1]]['cost'] for e in edge_list]
    cost_min = min(cost_list)
    cost_argminlist = [i for i,cost in enumerate(cost_list) if cost == cost_min]
    cost_argmin = random.choice(cost_argminlist)
    return edge_list[cost_argmin]

def get_scc_list(G):
    scclist = [scc for scc in nx.strongly_connected_components(G)
               if len(scc)>1]
    #scclist.reverse()
    return scclist

def remove_edges_in_scc(G,scc,edgecostfcn=rss_edgecost):
    """
    find a cycle in the strongly-connected-component (scc) of graph G
    remove the the lowest-cost edge in that cycle
    repeat, until a pair of edges fails to produce a cycle
    return the list of edges that were removed
    """
    elist=[]
    Gsub = nx.subgraph(G,scc).copy()  ## Gsub is frozen by default?
    while True: ## len(elist) < 50: ## True:
        #(a,b) = random.sample(scc,2)
        v.vprint("scc:",type(scc),len(scc))
        a, = random.sample(list(scc),1)
        #i = random.randint(0,len(scc)-2)
        #(a,b) = scc[i:i+2]
        try:
            #path01 = nx.shortest_path(G,source=a,target=b)  ## a->b
            #path10 = nx.shortest_path(G,source=b,target=a)  ## b->a
            #nodelist = path01[:-1] + path10[:-1]             ## a->b.pre + b->a.pre
            nodelist = my_bidirectional_shortest_path(Gsub,a)
            nodelist = nodelist[:-1]
            #print "a=",a,"nodelist=",len(nodelist),"elist=",len(elist)
            if not nodelist:
                raise RuntimeError("Empty path returned by shortest_path")
        except nx.NetworkXNoPath:
            if not elist:
                raise RuntimeError("No cycles found in scc -- "
                                   "should never happen!")
            return elist

        ecycle = list(zip(nodelist,nodelist[1:]+nodelist[:1])) ## list of edges
        e = mincost_edge_in_cycle(G,ecycle,edgecostfcn=edgecostfcn)
        elist.append(e)
        Gsub.remove_edge(*e)
        G.remove_edge(*e)
        v.vprint("Remove edge: ",e,
                 "from cycle of length",len(nodelist),
                 "in subgraph with",len(scc),"nodes")
    return elist


def decycle(G,str_edgecost='rss'):
    """
    Non-recursive (ie, iterative) version of decycle;
    more closely matches what is in the paper's pseudocode
    """
    tso=time.process_time()

    edges_removed = [] ## list of edges to be removed from toplevel graph GG

    if not str_edgecost:
        str_edgecost='rss'
    try:
        edgecost_fcn = edgecost_fcn_from_string[str_edgecost]
    except KeyError:
        print("Invalid heuristic edge-cost: [%s]; use [rss]"%str_edgecost)
        edgecost_fcn = rss_edgecost

    set_edgecost_all(G,edgecostfcn=edgecost_fcn)

    cnt_getscc = 1
    scclist = get_scc_list(G)
    while scclist:
        for scc in scclist:
            elist = remove_edges_in_scc(G,scc,edgecostfcn=edgecost_fcn)
            v.vprint("M:",cnt_getscc,"rm",len(elist),
                     "edges from scc with",
                     len(scc),"nodes")
            if len(elist)==0:
                print("scc=",scc)
                print("G=",nx.subgraph(G,scc))
                raise RuntimeError
            edges_removed.extend(elist)
        cnt_getscc += 1
        scclist = get_scc_list(G)

    print("DECYCLE: ",str_edgecost,"Removed",len(edges_removed),"edges",end="")
    totalcost=0
    maxcost=0
    for e in edges_removed:
        try:
            cost = G.nodes[e[0]]['count'] + G.nodes[e[1]]['count']
        except KeyError as err:
            v.print("e=",e,"; e[0],e[1]=",e[0],e[1])
            raise RuntimeError("KeyError") from err
        totalcost += cost
        maxcost = max([maxcost,cost])
    print("max=",maxcost,"total=",totalcost,end=" ")
    if len(edges_removed)>0:
        print("average=",totalcost/len(edges_removed))
    print("DECYCLE time:",time.process_time()-tso,"sec")

    return G
