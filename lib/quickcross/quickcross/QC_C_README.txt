Quick Cross heuristic for crossing minimisation in graphs.

An implementation of https://arxiv.org/abs/1804.09900 in the language C. 

Copyright (C) 2018 Michael Haythorpe, Alex Newcombe.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


COMPILING:

No packages are required except for the maths library in Linux.

Using GCC in Windows:

gcc QC.c -o QC.exe

Using GCC in Linux:

gcc -lm QC.c -o QC



RUNNING:

Run using a command line according to the options and examples listed below. Alternatively enter the command QC to see more detailed information on-screen.

Graphs can be input as either edge lists, or graph6 format. Alternatively we provide several built-in families of graphs which can be constructed inside the algorithm. Type the command QC to see all the options for these. Here is an example of constructing and solving a complete graph on 20 vertices:

QC -c 1 20

There are two options for input style and some options are only available for one of the styles.  First, if the graph is in graph6 format, you can input a graph using only the command line, for example:

QC -g Ys\_gC@?O@_E?B??_?G?A??K??o??W??C??@???O??@???C???F???F?

Second, you can point the algorithm to a file which has the graph inside as either graph6 format or a traditional edge list (coxeter_edgelist.txt is an example). For an edgelist format use the -f option, e.g:

QC -f coxeter_edgelist.txt

If the graph is in graph6 format, a file may have many graphs inside (one on each line) and the algorithm will solve and output all of them:

QC -g coxeter_graph6.txt

Initial embedding scheme options (-e):

1 : Kamada Kawai spring model drawing.
2 : embed onto circle.
3 : planar build embedding - build upon a chordless cycle, one vertex at a time.
4 : specify a crossing order list corresponding to a planarisation of a drawing of the graph (see QC_EXTRA_INFO.txt for details on this option).
5 : specify coordinates in the plane corresponding to vertex locations. The edges will then be drawn as straight lines and the graph will be solved from there (see QC_EXTRA_INFO.txt for details on this option).

Minimising scheme options (-m):

1 : use bigface scheme.
2 : use improvements as soon as they are found.
3 : search all vertices for best improvement before proceeding.

Output scheme (-o):

1 : only print the final number of crossings.
2 : print final number of crossings and a crossing order list. The crossings along an edge e=(u,v) where u<v, are listed in the order starting from v.
3 : print the edge list of a planar graph corresponding to the planarisation of the final solution. The original vertices are labelled 1,..,n and the planarised dummy vertices are labelled n+1,...,n+k.

Random seed (-r):

A natural number specifying the random stream. If set to zero then no randomisation occurs.

Stop when a good drawing is found (-s):

Terminates the run early if an embedding is found with no more than the specified number of crossings.

Upper bound (-u):

Continues to run the graph using different random seeds until an embedding with no more crossings than the upper bound is found. A random seed must be supplied here and serves as a starting point for the permutations. The random seed is incremented after each failure.

Verbosity (-v):

Style of the on-the-fly information during the iterations.

0 : only display final number of crossings.
1 : display number of crossings at each iteration.
2 : display all information.

Output options (-w):

If a filename is specified, the output is written to that file, otherwise output is printed on-screen.







