// Code was used from these sources:
// https://www.geeksforgeeks.org/greedy-algorithms-set-2-kruskals-minimum-spanning-tree-mst/
// https://www.topcoder.com/community/data-science/data-science-tutorials/range-minimum-query-and-lowest-common-ancestor/

#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <cmath>
#include <limits.h>

using namespace std;


typedef struct cartesian_node 
{
	// Struct that contains info for a Cartesian Node
	int id;
	int value;
	int parent; 
	int left=0;
	int right=0;
} Node;

typedef struct edge_struct 
{
	int u, v;
	int w;
} edge;

bool compare_edge(edge e1, edge e2) { return e1.w < e2.w; }

// To represent Disjoint Sets
struct DisjointSets
{
    int *parent, *rnk;
    int n;
 
    // Constructor.
	DisjointSets() { this->n = 0; }

    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];
 
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;
 
            //every element is parent of itself
            parent[i] = i;
        }
    }
 
    // Find the parent of a node 'u'
    // Path Compression
    int find(int u)
    {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }
 
    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);
 
        /* Make tree with smaller height
           a subtree of the other tree  */
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else // If rnk[x] <= rnk[y]
            parent[x] = y;
 
        if (rnk[x] == rnk[y])
            rnk[y]++;
    }
};

class Graph
{
	vector< vector<int> > adj_list;

public:
	int V, E;
	Graph(int V, int E);
	Graph() {}
	vector<edge> edges;

	void printArray(vector<vector<int>>);
	void addEdge(int u, int v, int d);
	void initGraph(int NumEdges);
	void kruskalMST(vector<edge> edges_list);
	void initQueries();
	vector<int> discover(int node, int came_from);
};

Graph::Graph(int V,int E)
{
	this->V = V;
	this->E = E;

	adj_list = vector< vector<int>>(V+1, vector<int>());
}

void Graph::addEdge(int u, int v, int d)
{
	adj_list[u].push_back(v);
	adj_list[v].push_back(u);

	edges.push_back({u, v, d});
}

void Graph::initGraph(int NumEdges)
{
	int i,j;
	int d;

	for (int k=0; k<NumEdges; k++) {
		cin >> i;
		cin >> j;
		cin >> d;

		// Put the dist d on pair(i,j)
		addEdge(i,j,d);
	}
}

void Graph::kruskalMST(vector<edge> edges_list)
{
    // Sort edges in increasing order on basis of cost
    sort(edges_list.begin(), edges_list.end(), compare_edge);
	
    // Create disjoint sets
    DisjointSets ds(V);
 
    // Iterate through all sorted edges
    vector<edge>::iterator it;
    for (it=edges_list.begin(); it!=edges_list.end(); it++)
    {
        int u = it->u;
        int v = it->v;
 		
        int set_u = ds.find(u);
        int set_v = ds.find(v);
	 
        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
            // Current edge will be in the MST
            // so print it
            addEdge(u, v, it->w); 
 
            // Merge two sets
            ds.merge(set_u, set_v);
        }
    }
}

class ParentInfo
{
	// This Class is an extra layer on top
	// of Union-Find Structure, since we 
	// want to make sure that each time we merge two
	// Nodes, we can know their parent
	public:
	vector<int> parent;
	DisjointSets ds;
		
	ParentInfo(int V) {
		parent = vector<int>(V+1);
		for (int i=1; i<=V; i++) { parent[i] = i; }

		DisjointSets ds(V);
		this->ds = ds;
	}

	void set_parent(int u, int v, int p)
	{
		// When merging make sure the new parent if both
		// nodes will be p
		ds.merge(u, v);
		parent[ds.find(u)] = p;
	}

	int get(int u)
	{
		// The parent of the Node will be the
		// parent of the parent Node of u
		return parent[ds.find(u)];
	}
};

class QueryGuru
{
	int V,E;
	Graph mst;
	Node* cartesian;
	int *Ev, *L, *H;
	int** M;

public:
	QueryGuru(Graph mst);
	void createCartesian();
	void createELH(int node, int level, int* counter);
	void preprocessQueries(int* A, int N);
	int RMQ(int* L, int i, int j);
	int query(int i, int j);
};

QueryGuru::QueryGuru(Graph mst)
{
	this->mst = mst;
	this->V = mst.V;
	this->E = mst.E;

	this->cartesian = (Node*)malloc((V + E + 1) * sizeof(Node));

	// Create the Cartesian tree from the mst
	createCartesian();

	// Prepare the arrays for the Euler Tour
	this->Ev = (int*)malloc(2*(V+E)*sizeof(int));
	this->L = (int*)malloc(2*(V+E)*sizeof(int));
	this->H = (int*)malloc(2*(V+E)*sizeof(int));

	for (int i=0; i<2*(V+E); i++) {
		Ev[i] = L[i] = H[i] = 0;
	}

	// Euler Tour to convert our problem to an RMQ
	int counter = 0;
	createELH(V + E, 0, &counter);

	int nV = 2*(V+E);
	this->M = (int**)malloc(nV*sizeof(int*));
	for (int i=0; i<nV; i++) { M[i] = (int*)malloc(ceil(log2(nV))*sizeof(int)); }

	// Sparce Matrix implementation of RMQ
	preprocessQueries(L, 2*(E+V));
}

void QueryGuru::createCartesian()
{

	//Sort the edges in increasing order
	sort(mst.edges.begin(), mst.edges.end(), compare_edge);

	// Initialize the vector containing the nodes
	//vector<Node> cartesian(V + E + 1);
	for (int i=1; i<=V+E; i++) {
		cartesian[i].id = i;
		cartesian[i].value = i;
		cartesian[i].left = 0;
		cartesian[i].right = 0;
		cartesian[i].parent = 0;
		// cartesian[i] = {i, i, 0, 0, 0};
	}

	// Create the struct that remembers the parents
	ParentInfo parents(V);

	for (int i=V+1; i<=V + E; i++) {
		edge min_edge = mst.edges[i-V-1];
		int u = min_edge.u, v = min_edge.v;

		// Put the correct children to 
		cartesian[i].left = parents.get(u);
		cartesian[i].right = parents.get(v);
		cartesian[i].value = min_edge.w;

		// Update the parent values in the cartesian tree
		cartesian[ cartesian[i].left ].parent = i;
		cartesian[ cartesian[i].right].parent = i;
		
		// Update the values of the parents for the nodes
		parents.set_parent(u, v, i);
	}
}

void QueryGuru::createELH(int node, int level, int* counter)
{
	// This function does the Eulerian Tour
	(*counter)++;

	// Update the values of the arrays
	Ev[*counter] = node;
	L[*counter] = level;
	if (H[node] ==  0) { H[node] = *counter; }

	// DFS to the left child
	if (cartesian[node].left != 0) {
		createELH(cartesian[node].left, level+1, counter);

		// Put the new entries
		(*counter)++;
		Ev[*counter] = node;
		L[*counter] = level;
	}

	if (cartesian[node].right != 0) {
		createELH(cartesian[node].right, level+1, counter);

		// Put the new entries
		(*counter)++;
		Ev[*counter] = node;
		L[*counter] = level;
	}
}

void QueryGuru::preprocessQueries(int* A, int N)
{
	// This is the Sparce Matrix Optimization Method
	// for solving RMQs with <O(nlogn), O(1)>
	// preprocess and query complexities
	int i, j;

	//initialize M for the intervals with length 1
	for (i = 0; i < N; i++)
		M[i][0] = i;
	//compute values from smaller to bigger intervals
	for (j = 1; 1 << j <= N; j++)
		for (i = 0; i + (1 << j) - 1 < N; i++)
			if (A[M[i][j - 1]] < A[M[i + (1 << (j - 1))][j - 1]])
				M[i][j] = M[i][j - 1];
			else
				M[i][j] = M[i + (1 << (j - 1))][j - 1];
}

int QueryGuru::RMQ(int* A, int i, int j)
{
	// This answer answers an RMQ in O(1) time
	// based on the M[n][logn]. This matrix needs
	
	int k = log2(j - i + 1);
	if (A[M[i][k]] <= A[M[j - (int)pow(2,k) + 1][k]]) 
		return M[i][k];
	else
		return M[j - (int)pow(2,k) + 1][k];
}

int QueryGuru::query(int i, int j)
{
	// This is the conversion of the LCA to RMQ problem, since
	// we can solve RMQ query in O(1), we can do the same with
	// the LCA on the cartesian tree
	return cartesian[ Ev[ RMQ(L, min(H[i],H[j]), max(H[i],H[j])) ] ].value;
}

int main()
{
	int V, E, Q, u, v, query;
	cin >> V;
	cin >> E;

	// init the Input Graph
	Graph g(V,E), mst(V, V-1);
	
	// Read from the Input
	g.initGraph(E);

	// Create the mst
	mst.kruskalMST(g.edges);
	
	QueryGuru Kimonas(mst);
	
	cin >> Q;
	for (int i=0; i<Q; i++) {
		cin >> u;
		cin >> v;
		query = Kimonas.query(u, v);
		cout << query << endl;
	}

	return 0;
}
