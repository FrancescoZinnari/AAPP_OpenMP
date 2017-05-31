#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BC.h"

int * nodes;
Edge * edges;
int N;
int M;



void buildGraph(char* fpathIn) {
	// read as adjacent list
    int u, v;
	double weight;

	FILE * f = fopen(fpathIn, "r");
    if(f == NULL){
        printf("ERROR: something went wrong when opening input file %s.",fpathIn);
        exit(1);
    }
	fscanf(f, "%d %d", &N, &M);
	int * adjNodes = (int *)malloc(N * sizeof(int));
	AdjEdge * adjEdges = (AdjEdge *)malloc(2 * M * sizeof(AdjEdge));
	memset(adjNodes, -1, N * sizeof(int));
	memset(adjEdges, 0, M * sizeof(AdjEdge));



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

/*
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

*/
}

void outputBC(double* arrayBC, char *fpathOut,int N){
    FILE *f = fopen(fpathOut,"w");
    if(f == NULL){
        printf("ERROR: something went wrong when opening output file %s.",fpathOut);
        exit(1);
    }

    for(int i=0; i<N;i++){
        fprintf(f,"%d,%.1f\n", i, arrayBC[i]);
    }
}


int main(int argc, char *argv[]) {
    double *arrayBC;

    if(argc<3){
        printf("USAGE: you need to specify (1) path to input graph file and (2) path to output BC file");
        exit(1);
    }

	buildGraph(argv[1]);
    arrayBC = BC(N, M, nodes, edges);
    outputBC(arrayBC,argv[2],N);


    for(int i=0;i<N;i++){
        printf("\nBC of node %d -> %f",i,arrayBC[i] );
    }

    free(arrayBC);
	return 0;
}
