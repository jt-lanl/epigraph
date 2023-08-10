## Modified from source code in NetworkX package

import networkx as nx
import heapq
import sys


def my_bidirectional_shortest_path(G,source,target=None):
    """
    Modified version of bidirectional_shortest_path that interprets
    target=None to mean find a cycle that includes the source=target node

    Return a list of nodes in a shortest path between source and target.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       starting node for path

    target : node label
       ending node for path

    Returns
    -------
    path: list
       List of nodes in a path from source to target.

    Raises
    ------
    NetworkXNoPath
       If no path exists between source and target.

    See Also
    --------
    shortest_path

    Notes
    -----
    This algorithm is used by shortest_path(G,source,target).
    """

    if source == target:
        return [source]

    if target is None:
        target=source

    # call helper to do the real work
    pred,succ,w=_mybidirectional_pred_succ(G,source,target)

    # build path from pred+w+succ
    path=[]

    # from w to target
    path.append(w)
    while w != target:
        w = succ[w]
        path.append(w)

    # from source to w
    w = pred[path[0]]
    path.insert(0,w)
    while w != source:
        w = pred[w]
        path.insert(0,w)

    return path

def _mybidirectional_pred_succ(G, source, target):
    """Bidirectional shortest path helper.

       Returns (pred,succ,w) where
       pred is a dictionary of predecessors from w to the source, and
       succ is a dictionary of successors from w to the target.
    """
    # does BFS from both source and target and meets in the middle
    # If target == source, we seek a cycle

    # handle either directed or undirected
    if G.is_directed():
        Gpred=G.predecessors
        Gsucc=G.successors
    else:
        Gpred=G.neighbors
        Gsucc=G.neighbors

    # predecesssor and successors in search
    succ={}
    pred={}

    # initialize fringes, start with forward
    forward_fringe=[source]
    reverse_fringe=[target]

    while forward_fringe and reverse_fringe:
        if len(forward_fringe) <= len(reverse_fringe):
            this_level=forward_fringe
            forward_fringe=[]
            for v in this_level:
                for w in Gsucc(v):
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w]=v
                    if w in succ or w==target:  
                        return pred,succ,w # found path
        else:
            this_level=reverse_fringe
            reverse_fringe=[]
            for v in this_level:
                for w in Gpred(v):
                    if w not in succ:
                        succ[w]=v
                        reverse_fringe.append(w)
                    if w in pred or w==source:  
                        return pred,succ,w # found path


    raise nx.NetworkXNoPath("No path between %s and %s." % (source, target))

