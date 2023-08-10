# EpiGraph: optimal epitope coverage

## Algorithm and Purpose

*EpiGraph* is a graph-based algorithm for optimizing epitope
coverage. By "epitope" we really mean "potential T-cell epitope" and
by that we mean a _k_-mer of amino acids, and we almost always take
_k=9_ (though _k=8_ to _k=12_ are also reasonable choices). Given a
set of natural sequences corresponding to the natural variability of a
protein in the wild, the goal is to design a single sequence, or a
cocktail of several sequences, such that the potential epitopes (that
is, all _k_-character sub-sequences) found in the designed sequences
cover as many as possible of the epitopes in the natural sequences.

The algorithm was introduced in [1], and is more carefully described
in [2], including a discussion of the various heuristics used for the
de-cycle step.  The *EpiGraph* algorithm shares goals with the
*Mosaic* [3] vaccine design algorithm, except that *Mosaic* uses a
probabilistic genetic algorithm to optimize coverage, while *EpiGraph*
uses a graph-based approach to maximize the same potential epitope
coverage that the Mosaic algorithm maximizes, but with better
computational efficiency; and, under some conditions, with provably
optimal solutions.

The graph that *EpiGraph* uses is a deBrujn graph in which each node
is a _k_-mer of amino acids. Two nodes are connected when the last
_k-1_ amino acids of one node agree with the first _k-1_ amino acids
of the other.  Each node is assigned a value according to the number
of sequences in the input dataset that contain the associated _k_-mer.
A path through this directed graph is scored according to the sum of
the values of all of the nodes in the path.  Since we do not want to
count a node's value more than once, we seek paths that do not
re-visit old nodes.  If the original path has no loops, then it is a
directed acyclic path (DAG), and for such a graph finding the path
with the highest value can be efficiently accomplished.

If there are no cycles in the deBrujn graph that is generated from the sequences, then *EpiGraph* produces a provably optimal solution.  In practice (and this obviously depends on the application):
* There are usually a few cycles are produced in the initial graph.
* The heuristic de-cycling is relatively efficient in removing those cycles.
* Since de-cycling has a random component, multiple de-cycling trials can be pursued to improve performance, though the difference between trials is usually very slight.
* *EpiGraph* is usually much much faster than *Mosaic* and usually
achieves higher coverage

## Web implementation

The *EpiGraph* algorithm can be run directly from the web (<https://www.hiv.lanl.gov/content/sequence/EPIGRAPH/epigraph.html>). The algorithm on the web is implemented using the software in this package.

## Python implementation

This *EpiGraph* software package is a python implementation of the epigraph algorithm. It takes as input a fasta file (several other formats are supported, as provided in the *sequtil* package) of amino acid sequences that characterizes the natural variability of the protein whose epitopes you want to cover.  It provides on output a fasta (or, possibly, some other format) file with sequences of artificial proteins whose epitopes cover (*ie*, match) as many of the epitopes as it can of the epitopes in the input sequences.

### Files

* `epiun.py`: main routine for running epigraph algorithm on unaligned sequences
* `epiuutil.py`: various utility routines called by `epiun.py`
* `episequtil.py`: various sequence-based utilities
* `decycle.py`: routines for removing cycles from the directed graph, and producing a directed acyclic graph (DAG)
* `config.py`: contains global variables
* `nxutil.py`: utilities for interacting with the networkx package
* `newcycle.py`: re-implementation of the shortest path algorithm

#### Auxiliary files

These files are also included in the package, but are generic (and are included in several other packages as well)

* `verbose.py`: routines for writing information message if a verbose flag is set
* `sequtil.py`: routines for reading and writing fasta files (several other sequence formats are also supported)

### Typical usage

    
#### Simple run

    python -m epiun -i input.fasta -E9 -M3 -out output.fasta
	
Here `input.fasta` is a list of sequences (they can be aligned or not, but the algorithm does not make use of the alignment, just of the sequences), `-E9` indicates that 9-mers will be used as the potential T-cell epitopes, `-M3` indicates that a cocktail of three sequences will be output into the file `output.fasta`.

#### The 'help' command-line option

First of all, note that the `--help` (or `-h`) command line option produces a list of available options; for example, if you type

    python epiun.py -h
	
then you'll get something like the following:

  
    usage: epiun.py [-h] [-E E] [--inseq INSEQ] [--out OUT]
                    [--outseqnames OUTSEQNAMES] [-T T] [--seed SEED]
                    [--Mxtra MXTRA] [-M M] [--weightsfile WEIGHTSFILE] [--tmr TMR]
                    [--reverse] [--nodelowcount NODELOWCOUNT]
                    [--edgelowcount EDGELOWCOUNT] [--rmweaknodes] [--dag DAG]
                    [--writegraph WRITEGRAPH] [--usegraph USEGRAPH]
                    [--fixinit FIXINIT] [--hintinit HINTINIT] [--verbose]
    
    epiun: EPIgraph UNaligned
    
    optional arguments:
      -h, --help            show this help message and exit
      -E E                  Epitope length
      --inseq INSEQ, -i INSEQ
                            input sequence filename
      --out OUT             output file with vaccine sequences
      --outseqnames OUTSEQNAMES
                            comma-separated list of names of vaccine sequences
      -T T                  Number of trials
      --seed SEED           random seed
      --Mxtra MXTRA         extra iterations beyond M
      -M M                  Number of seqs in cocktail
      --weightsfile WEIGHTSFILE
                            File with weights associated to names
      --tmr TMR             (too many repeats) remove sequences with more than TMR repeats
      --reverse             reverse sequences
      --nodelowcount NODELOWCOUNT
                            delete nodes with count <= LOWCOUNT
      --edgelowcount EDGELOWCOUNT
                            delete edges whose nodes have low counts
      --rmweaknodes         delete nodes whose self and all neighbors have count=1
      --dag DAG             Heuristic for dag-ifying graph: rss, sum, etc
      --writegraph WRITEGRAPH
                            write graph object to python pickle file
      --usegraph USEGRAPH   read graph object from python pickle file
      --fixinit FIXINIT     file with fixed initial vaccine sequences
      --hintinit HINTINIT   file with suggested initial vaccine sequences
      --verbose, -v         verbose

# NetworkX

*EpiGraph* makes extensive use of the *NetworkX* (<https://networkx.org/>) package for network algorithms, and we are grateful to the authors of that package for making it open-source. The *EpiGraph* package also includes (in the `newcycle.py` file) a re-implementation of the "bidirectional_shortest_path" function that enables it to find paths from a node to itself (if such paths exist; these are cycles, and this is one step in our ultimate aim of removing all cycles in the graph).

# References

[1] J. Theiler, H. Yoon, K. Yusim, L. J. Picker, K. Frueh, and B. Korber. "Epigraph: A vaccine design tool applied to an HIV therapeutic vaccine and a pan-Filovirus vaccine." Scientific Reports 6 (2016) 33987. 

[2] J. Theiler and B. Korber. "Graph-based optimization of epitope coverage for vaccine antigen design." Statistics in Medicine 37 (2017) 181-194. 

[3] W. Fischer, S. Perkins, J. Theiler, T. Bhattacharya, K. Yusim, R. Funkhouser, C. Kuiken, B. Haynes, N. L. Letvin, B. D. Walker, B. H. Hahn, and B. T. Korber. "Polyvalent vaccines for optimal coverage of potential T-cell epitopes in global HIV-1 variants." Nature Medicine 13 (2007) 100-106.

# COPYRIGHT ASSERTION: C23039-EpiGraph

(c) 2023. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

# LICENSE (GPLv3)

This software is open source under the GNU General Public License ; version 3.0 of the License. See file `GPLv3.pdf` for details.