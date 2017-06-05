#ifndef BC_INCLUDED
#define BC_INCLUDED


typedef struct {
	int v, next;
	double weight;
}AdjEdge;

typedef struct {
	int v;
	double weight;
} Edge;

typedef struct Neigh *Neighbours;
struct Neigh{
    Edge* neighbours;
    int nOfNeighbours;
};

double* BC(int N,int M,int* nodes,Edge* edges /*calNode, calEdge, bcpath.c_str(), ebcpath.c_str(), warp_size*/);

#endif // BC_INCLUDED
