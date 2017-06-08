#include "BC.h"
#include <omp.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

/*** Returns a Neighbours type (pointer to struct) containing neighbours of node "node")***/

Neighbours neighboursOfNode(int node,int N, int M, int* nodes, Edge* edges){
    int begin;
    int end;
    Neighbours x;

    begin = nodes[node];
    end = nodes[node+1];

    if(end-begin==0){
        return NULL;
    }

    x = malloc(sizeof(*x));
    x->neighbours = &edges[begin];
    x->nOfNeighbours = end-begin;

    return x;
}


/*** Returns minimum edge weight from node "node"***/

double computeMinEdge(int node, int N, int M, int* nodes, Edge* edges){
    Neighbours nn = neighboursOfNode(node,N,M,nodes,edges);
    double minEdgeWeight= DBL_MAX;

    if(nn==NULL || nn->nOfNeighbours <= 0){
        fprintf(stderr,"Warning: node %d does not have neighbours\n", node);
        return DBL_MAX;
    }

    for(int i=0; i< (nn->nOfNeighbours);i++){
        if(minEdgeWeight > nn->neighbours[i].weight){
            minEdgeWeight = nn->neighbours[i].weight;
        }
    }

    return minEdgeWeight;

}


/*** Computes dependencies from single node "root" and update BC***/

void processSingleRoot(int root, int N, int M, int* nodes, Edge* edges, double* arrayBC, double *minEdge){
    int *unsettled, *frontier, *ends, *settled;
    double *distance, *nMinPath, *dep;
    double delta;
    int flen,slen,endsLen, depth, start, end;

    unsettled = malloc(N*sizeof(int));
    frontier = malloc(N*sizeof(int));
    ends = malloc(N*sizeof(int));
    settled = malloc(N*sizeof(int));
    distance = malloc(N*sizeof(double));
    nMinPath = malloc(N*sizeof(double));
    dep = malloc(N*sizeof(double));


    {
        /*** 1. INITIALIZATION ***/

        for(int i=0; i<N; i++){
            unsettled[i] = 1;
            distance[i] = DBL_MAX;
            nMinPath[i] = 0.0;
            dep[i] = 0.0;
            ends[i] = 0;
            settled[i] = 0;
        }

        distance[root] = 0;
        nMinPath[root] = 1;
        unsettled[root] = 0;
        frontier[0] = root;
        flen = 1;
        settled[0] = root;
        slen = 1;
        ends[0] = 0;
        ends[1] = 1;
        endsLen = 2;
        delta = 0;

        /*** 2. DIJKSTRA ***/

        while(delta < DBL_MAX)
        {
            /*** TENTATIVE DISTANCE ***/

            for(int i=0; i<flen;i++){
                int x = frontier[i]; //<--- node from frontier

                //acquire neighbours of x
                Neighbours nx = neighboursOfNode(x, N, M, nodes, edges);

                if(nx==NULL || nx->nOfNeighbours <= 0){
                    fprintf(stderr,"Unexpected behaviour, node %d should have neighbours", x);
                    exit(1);
                }

                //for each neighbour
                for(int j=0; j<(nx->nOfNeighbours); j++){
                    Edge ey = nx->neighbours[j]; //<--- edge to neighbour node
                    int y = ey.v;

                    if (unsettled[y]==1 && (distance[y] > (distance[x] + ey.weight))){
                        distance[y] = distance[x] + ey.weight;
                        nMinPath[y] = 0;
                    }

                    if(distance[y] == distance[x] + ey.weight){
                        nMinPath[y] = nMinPath[y] + nMinPath[x];
                    }

                } //END FOR ON NEIGHBOURS
            } //END FOR ON FRONTIER



            /*** COMPUTE DELTA ***/

            delta = DBL_MAX;

            for(int i=0; i<N;i++){
                if(unsettled[i] && distance[i]<DBL_MAX){
                    if(distance[i] + minEdge[i] < delta){
                        delta = distance[i] + minEdge[i];
                    }
                }
            }

            /*** UPDATE FRONTIER ***/

            flen = 0;

            for (int i=0; i<N; i++){
                if(unsettled[i] && distance[i]<delta){
                    unsettled[i]=0;
                    int t;
                    t = flen;
                    flen = flen+1;
                    frontier[t]=i;

                }
            }


            /*** UPDATE DAG ***/

            if(flen>0){

                ends[endsLen] = ends[endsLen - 1] + flen;
                endsLen++;

                for (int i=0; i<flen; i++){
                    settled[slen+i]=frontier[i];
                }
                slen = slen + flen;
            }


        } //END OF DIJKSTRA WHILE LOOP


        /*** 3. DEPENDENCY SUMMATION ***/

        depth = endsLen - 1;
        while(depth > 0){
            start = ends[depth -1];
            end = ends[depth] ;

            for(int i=0; i < end - start; i++){
                int w = settled[start+i];

                if(w==root)
                    continue;

                double dsw = 0;
                Neighbours nw = neighboursOfNode(w,N,M,nodes,edges);

                if(nw==NULL || nw->nOfNeighbours <= 0){
                    fprintf(stderr,"Unexpected behaviour, node %d should have neighbours", w);
                    exit(1);
                }

                for(int j=0; j < (nw->nOfNeighbours); j++){
                    Edge ev = nw->neighbours[j];
                    int v = ev.v;

                    if(distance[v] == distance[w]+ev.weight){
                        dsw = dsw + (nMinPath[w]/nMinPath[v])*(1+dep[v]);
                    }
                }

                dep[w] = dsw;

                if (w!=root){
                    #pragma omp atomic
                        arrayBC[w] = arrayBC[w]+dep[w];
                }
            }
            depth = depth - 1;
        } //end of while on DAG


    } //END OF PARALLEL REGION

    free(unsettled);
    free(frontier);
    free(ends);
    free(settled);
    free(distance);
    free(nMinPath);
    free(dep);
}


double* BC(int N,int M,int* nodes,Edge* edges){
    double* arrayBC, *minEdge;
    arrayBC = malloc(N*sizeof(double));
    minEdge = malloc(N*sizeof(double));


    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0;i<N;i++){
            arrayBC[i] = 0;
            minEdge[i] = computeMinEdge(i,N,M,nodes,edges);
        }

        #pragma omp for
        for(int i=0; i<N; i++){
            if(minEdge[i]!=DBL_MAX){
                processSingleRoot(i,N,M,nodes,edges,arrayBC,minEdge);
            }
        }
    }

    free(minEdge);
    return arrayBC;
}

