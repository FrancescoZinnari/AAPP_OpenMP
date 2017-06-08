#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

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
}

void outputBC(double* arrayBC, char *fpathOut){
    FILE *f = fopen(fpathOut,"w");
    if(f == NULL){
        printf("ERROR: something went wrong when opening output file %s.",fpathOut);
        exit(1);
    }

    for(int i=0; i<N;i++){
        fprintf(f,"%d,%f\n", i, arrayBC[i]);
        //printf("\nBC of node %d -> %f",i,arrayBC[i] );
    }

    fclose(f);
}

int main(int argc, char *argv[]) {
    double *arrayBC;
    char* inputPath, *outputPath;
    unsigned int num_threads;
    double start, end;

    if(argc<4){
        printf("USAGE: you need to specify (1) path to input graph file, (2) path to output BC file and (3) number of threads of the pool");
        exit(1);
    }

    inputPath = argv[1];
    outputPath = argv[2];
    num_threads = (unsigned int) atoi(argv[3]);

    #ifdef _OPENMP
        /* Set the number of threads */
        omp_set_num_threads(num_threads);
    #endif

	buildGraph(inputPath);
	printf("Finished acquiring graph, now computing BC\n");

    start = omp_get_wtime();
    arrayBC = BC(N, M, nodes, edges);
    end =  omp_get_wtime();
    outputBC(arrayBC,outputPath);
    printf("Time elapsed %f", (end-start));

    free(arrayBC);
	return 0;
}
