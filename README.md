# MST-Cartesian-Trees
-
The other day in my CS Classes they gave us an assignment where given a graph with edges
and weights, we had to answer a range of queries for the Minimax Path Problem. The number of
queries was ~ 100.000 and the Graph could have up to V=30.000 Nodes. This means that a naive
implementation of just finding an MST and then with a modified DFS to calculate the minimax
distance between all the pairs wouldn't do the job since it would have O(V^2) complexity.

# On to Cartesian Trees
-
The most efficient solution for the problem would be to create a Cartesian Tree based from
the Minumum Spanning Tree which allows us to find the minimax distance between any pair of 
vertices to be queried in constant time per query, using lowest common ancestor queries in 
a Cartesian tree.

The only problem was that both Wikipedia and and some papers I found olnline with a basic 
Google Search described an Algorithm for creating the Cartesian with O(n) complexity but
by using a Structure for keeping track of decremental components of a tree. For this struct
there was of course no code available and only one paper describing it. So I had to tackle
the creation of the Cartesian with a different way, given that there were just two days remaining
for the assignment. 

In the end, I came up with a way to create the Cartesian Tree from the MST with an O(nlogn) 
complexity, which was pretty acceptable given that I had to do Kruskal in order to find the
actual MST. So in this repo I give the code of the full assignment which includes both the 
code for creating the Cartesian Tree, but also the code for transforming the LCA to RMQ for 
actually answering the queries.

Bellow I'll give a basic explanation of the algorithm used to create the Cartesian Tree.

# Creating the Cartesian Tree
-
To do this instead of creating the tree top-down we will create it Bottom-Up. We will be 
constantly changing the parents of the current nodes and keep creating the tree until we
get to the top. The top Node will be the node V + E (where E = V - 1) and it's value will
be the weight of the heaviest edge and the leaves will be the Nodes of the initial Graph. 
The resulting tree will be stored in an array of length V+E where each element is a type
CartesianNode Class which keeps track of its parent Node, as well as its children.

Firstly we will need a struct that can keep track of the parents of the nodes that we will
be changing. For this we will use a Union-Find Struct (hence the O(nlogn) complexity) but 
will a minor modification for keeping track of the parent Nodes.


With this, we will loop for each edge in sorted edges (with increasing order based on their 
weights) of the MST. In each iteration we will add another edge as a Node in the tree and
update the values of that node (its children) as well as its children (change the value of
their parent). 

# Answering the Queries
-
I won't get into the technical details here, since there is this [tons of helpful](https://www.topcoder.com/community/data-science/data-science-tutorials/range-minimum-query-and-lowest-common-ancestor/)
explanation of the procedure from here on. The general idea is that from this Cartesian 
tree we can find the minimax shortest path by finding the Lowest Common Ancestor in the 
Cartesian Tree for the two Nodes. But, finding an LCA is equivalent with an RMQ and since
by doing a O(nlogn) preprocess we can then answer in O(1) any RMQ, we can finaly also 
answer any query of the initial problem in constant time.
