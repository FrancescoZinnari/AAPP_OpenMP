#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BC.h"

int * nodes;
Edge * edges;
int N;
int M;



void buildGraph(char* fpath) {
	// read as adjacent list
	FILE * f = fopen(fpath, "r");
	fscanf(f, "%d %d", &N, &M);
	int * adjNodes = (int *)malloc(N * sizeof(int));
	AdjEdge * adjEdges = (AdjEdge *)malloc(2 * M * sizeof(AdjEdge));
	memset(adjNodes, -1, N * sizeof(int));
	memset(adjEdges, 0, M * sizeof(AdjEdge));
	int u, v;
	double weight;
	for (int i = 0; i < M; i++) {
		fscanf(f, "%d %d %lf", &u, &v, &weight);
		adjEdges[i].weight = weight;
		adjEdges[i].v = v;
		adjEdges[i].next = adjNodes[u];
		adjNodes[u] = i;
		adjEdges[i + M].weight = weight;
		adjEdges[i + M].v = u;
		adjEdges[i + M].next = adjNodes[v];
		adjNodes[v] = i + M;
	}
	fclose(f);

	//convert adjacent list to CRS format
	nodes = (int *)malloc((N+1) * sizeof(int));
	edges = (Edge *)malloc(2 * M * sizeof(Edge));
	int edgeLen = 0;
	nodes[0] = 0;
	for(int i = 0; i < N; i++) {
		int edgeIdx = adjNodes[i];
		while(edgeIdx != -1) {
			AdjEdge adjEdge = adjEdges[edgeIdx];
			edges[edgeLen].v = adjEdge.v;
			edges[edgeLen].weight = adjEdge.weight;
			edgeLen++;
			edgeIdx = adjEdges[edgeIdx].next;
		}
		nodes[i+1] = edgeLen;
	}
	free(adjNodes);
	free(adjEdges);


	//print nodes - DEBUG
	printf("NODES (COUNTER)\n");
	for(int i=0;i<=N;i++){
        printf("%d ---> %i\n",i,nodes[i]);
	}

    //print edges - DEBUG
	printf("EDGES\n");
	for(int i=0;i<2*M;i++){
        printf("%d ---> %d , %f\n",i,edges[i].v,edges[i].weight);
	}


}



int main(int argc, char *argv[]) {
    char *fpath = "/home/fzinnari/CB_workspace/OpenMP_AAPP/input/graph-simple.txt";
    double *arrayBC;

	buildGraph(fpath);
    arrayBC = BC(N, M, nodes, edges);

    for(int i=0;i<N;i++){
        printf("\nBC of node %d -> %f",i,arrayBC[i] );
    }

	return 0;
}
