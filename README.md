# Minimum Spanning Tree Interdiction Solver
This repo contains the code and problem instances accompanying the (forthcoming) paper "Interdiction of minimum spanning trees and other matroid bases" by Noah Weninger and Ricardo Fukasawa.

The solver is a standalone program which accepts minimum spanning tree (MST) interdiction or minimum MST edge blocker instances in a plain text format.
It then solves the instance to optimality using a combinatorial branch-and-bound algorithm and outputs the solution in a plain text format.

## Building
The code has only been tested on Linux (no support will be provided for other platforms).
There are no external dependencies besides a C++ compiler (new enough to support C++17) and CMake for building.
Although GCC or Clang can be used for building, it is strongly recommended to use Clang because it does a significantly
better job at optimizing the code: often, the speedup when using Clang is around 4x - 8x.
To get started, clone this repository or download the code. From the root directory of the repo, do the following:
```
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
$ make
```
This should produce an executable called `mstisolver`. If you have difficulties building or using the solver,
feel free to contact Noah (contact details can be found in the paper).

## Usage
The solver accepts instances in a simple text format, which is described in the **Input Format** section below.
We assume that you have extracted `instance.tar.zst` to a directory called `instances`, by using the following command (or similar)
in the root directory of the repo:
```
$ tar --zstd -xvf instances.tar.zst
```
Many instances are provided in the `instances` directory, with the file extension `.msti`. Files with the `.ans` extension
contain a single integer describing the optimal solution weight for the corresponding `.msti` file. If the `.ans` file
is missing, it means that the instance could not be solved within 1 hour by the test machine used in the paper.

For example, to solve the instance `instances/small/10.5.1.msti`:
```
$ build/mstisolver -f instances/small/10.5.1.msti
```
The output produced is
```
max_spanning_tree 7042
min_cut 757
canonicalize_time 0.050866
min_spanning_tree 1378
greedy_time 0.161016
greedy_lower_bound 2729
sol_val 2729
X 0000001100111000100000
Y 1111010001000000010101
status ok
wallclock_time 1.03077
cpu_time 0.977
num_bnb_nodes 202
num_bnb_pruned 0
```
It can be verified that 2729, the number following `sol_val`, matches the corresponding `.ans` file.
A full description of the output fields is described in the **Output Format** section below.
A number of command line options are available to configure various parameters of the solver.
For the full list, run
```
$ build/mstisolver --help
```
There is also a simple testing script included called run_tests.py, which can be used to solve many instances and
record timing information into a JSON file. To use the script, you can edit the code at the bottom
of the file to use the desired solvers, tests, and JSON file.

## Input Format

Instances are described by a simple plain text format.
The solver supports two types of instances: MST interdiction (`msti`) and minimum edge blocker (`mebsp`).
The first line of the instance file should contain
```
problem_type XXX
```
where XXX is either `msti` or `mebsp`.
After this, the number of vertices, number of edges, and C (for `msti` instances) or R (for `mebsp` instances) value should be specified.
The vertices and edges are specified by two lines:
```
n_verts X
n_edges Y
```
where X is the number of vertices and Y is the number of edges.
For `msti` instances, the capacity value C can be specified as follows:
```
cap C
```
For `mebsp` instances, the target weight value R can be specified as follows:
```
target_weight R
```
It is not strictly necessary to specifiy C or R in the instance file, as it can be overridden on the command line when invoking the solver.
Finally, the file should contain a number of lines specifying the edges of the graph. Each line has the form
```
edge s t w c
```
which indicates that there is an edge between vertices s and t with weight w and cost c.
The total number of such lines must equal the value of `n_edges`.
The edges must be specified in order of non-decreasing weight w.

## Output Format

By default, the solver outputs some (self-explanatory) extra information to standard error, which can be disabled with the `-q` flag.
In this section, we assume the `-q` flag is specified.
Below, we show the output from solving an instance, where each line of output is annotated with the meaning:
```
$ build/mstisolver -f instances/new/complete-n20-d1-cf1-c1.w10000.msti -q
max_spanning_tree 176726 # the maximum spanning tree weight of the graph
min_cut 19 # the global minimum cut cost of the graph
canonicalize_time 0.182418 # the time required for instance reduction using a minimum cost cocircuit
min_spanning_tree 9386 # the minimum spanning tree weight of the graph
dp_upper_bound_0 43256 # the weight of the DP upper bound with p=0
dp_upper_bound_time_0 7.39384 # the time in milliseconds needed to compute the DP upper bound with p=0
dp_upper_bound_1 43084 # the weight of the DP upper bound with p=1
dp_upper_bound_time_1 12.5034 # the time in milliseconds needed to compute the DP upper bound with p=1
greedy_time 31.5881 # the time in milliseconds needed to compute the greedy lower bound
greedy_lower_bound 31328 # the weight of the greedy lower bound
dp_upper_bound_2 42180 # the weight of the DP upper bound with p=2
dp_upper_bound_time_2 20.5235 # ...
dp_upper_bound_3 40877
dp_upper_bound_time_3 25.7584
dp_upper_bound_4 40681
dp_upper_bound_time_4 38.4542
dp_upper_bound_5 40066
dp_upper_bound_time_5 75.2055
dp_upper_bound_6 39394
dp_upper_bound_time_6 146.776
dp_upper_bound_7 39394
dp_upper_bound_time_7 285.352
dp_upper_bound_8 39082
dp_upper_bound_time_8 564.307
dp_upper_bound_9 38648
dp_upper_bound_time_9 1136.32
sol_val 32412 # the weight of the optimal solution
X 1111101111110111101100000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 # the optimal solution vector for the leader
Y 0000010000001000010011011000111110101000100001100000000000000000000001000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 # the optimal solution vector for the follower
status ok # indicates that there was no error encountered
wallclock_time 3512.87 # the total wall-clock time in milliseconds needed to solve the instance
cpu_time 7014.9 # the total CPU time in milliseconds needed to solve the instance
num_bnb_nodes 30440835 # the number of branch and bound nodes evaluated
num_bnb_pruned 14351804 # the number of times a node was pruned in branch and bound (not including children of pruned nodes)
```
