#include "BC.h"
#include <omp.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

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


double computeMinEdge(int node, int N, int M, int* nodes, Edge* edges){
    Neighbours nn = neighboursOfNode(node,N,M,nodes,edges);
    double minEdgeWeight= DBL_MAX;

    if(nn==NULL || nn->nOfNeighbours <= 0){
        fprintf(stderr,"Unexpected behaviour, node %d should have neighbours", node);
        exit(1);
    }

    //#pragma omp for reduction (min:minEdgeWeight)
    for(int i=0; i< (nn->nOfNeighbours);i++){
        if(minEdgeWeight > nn->neighbours[i].weight){
            minEdgeWeight = nn->neighbours[i].weight;
        }
    }


    return minEdgeWeight;

}

void processSingleRoot(int root, int N, int M, int* nodes, Edge* edges, double* arrayBC){
    int *unsettled, *frontier, *ends, *settled;
    double *distance, *nMinPath, *dep, *minEdge;
    double delta;
    int flen,slen,endsLen, depth, start, end;
    omp_lock_t* lock;

    unsettled = malloc(N*sizeof(int));
    frontier = malloc(N*sizeof(int));
    ends = malloc(N*sizeof(int));
    settled = malloc(N*sizeof(int));
    distance = malloc(N*sizeof(double));
    nMinPath = malloc(N*sizeof(double));
    dep = malloc(N*sizeof(double));
    minEdge = malloc(N*sizeof(double));
    lock = malloc(N*sizeof(omp_lock_t));



    #pragma omp parallel
    {
        /*** 1. INITIALIZATION ***/

        #pragma omp for
        for(int i=0; i<N; i++){
            unsettled[i] = 1;
            //frontier[i] = 0;
            distance[i] = DBL_MAX;
            nMinPath[i] = 0.0;
            dep[i] = 0.0;
            ends[i] = 0;
            settled[i] = 0;
            omp_init_lock(&(lock[i]));
            minEdge[i] = computeMinEdge(i,N,M,nodes,edges);
        }

        #pragma omp single
        {
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
        }

        /*** 2. DIJKSTRA ***/

        while(delta < DBL_MAX)
        {
            /*** TENTATIVE DISTANCE ***/

            #pragma omp for
            for(int i=0; i<flen;i++){
                //printf("Thread for index %d in for\n",i);
                int x = frontier[i]; //<--- node from frontier
                //printf("Thread for index %d, node %d\n", i, x);
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

                    //printf("neighbour of %d -> %d\n", x,y);

                    //acquire lock on y
                    omp_set_lock(&(lock[y]));
                    /*
                    printf("distance of x: %f\n",distance[x]);
                    printf("distance of y: %f\n",distance[y]);
                    printf("distance between x-y: %f\n",ey.weight);
                    printf("unsettled of y: %d\n",unsettled[y]);
                    */

                    if (unsettled[y]==1 && (distance[y] > (distance[x] + ey.weight))){
                        distance[y] = distance[x] + ey.weight;
                        nMinPath[y] = 0;
                        //printf("new distance for node %d: %f\n", y,distance[y]);
                    }

                    if(distance[y] == distance[x] + ey.weight){
                        nMinPath[y] = nMinPath[y] + nMinPath[x];
                    }

                    //release lock on w
                    omp_unset_lock(&(lock[y]));

                } //END FOR ON NEIGHBOURS
            } //END FOR ON FRONTIER



            /*** COMPUTE DELTA ***/

            #pragma omp single
            {   /*
                for(int i=0;i<N;i++){
                    printf("node %d, distance %f, nMinPath %f\n", i,distance[i], nMinPath[i]);
                }
                for(int i=0;i<N;i++){
                    printf("node %d, minEdge %f\n", i,minEdge[i]);
                }*/

                delta = DBL_MAX;
            }

            #pragma omp for reduction(min:delta)
            for(int i=0; i<N;i++){
                if(unsettled[i] && distance[i]<DBL_MAX){
                    if(distance[i] + minEdge[i] < delta){
                        delta = distance[i] + minEdge[i];
                        //printf("Delta: node %d having distance %f and minEdge %f\n",i,distance[i],minEdge[i]);
                    }
                }
            }

            /*** UPDATE FRONTIER ***/

            #pragma omp single

            {
                //printf("DELTA = %f\n",delta);
                flen = 0;
            }

            #pragma omp for
            for (int i=0; i<N; i++){
                if(unsettled[i] && distance[i]<DBL_MAX){
                    unsettled[i]=0;
                    int t;
                    #pragma omp critical
                    {
                        t = flen;
                        flen = flen+1;
                    }
                    frontier[t]=i;
                }
            }


            /*** UPDATE DAG ***/

            if(flen>0){
                 #pragma omp single nowait
                {
                    for(int i=0;i<flen;i++){
                        //printf("node of frontier pos %d ---> %d\n",i,frontier[i]);
                    }

                    ends[endsLen] = ends[endsLen - 1] + flen;
                    endsLen++;
                }
                #pragma omp for
                for (int i=0; i<flen; i++){
                    settled[slen+i]=frontier[i];
                }

                #pragma omp single
                {
                    slen = slen + flen;
                }
            }

            #pragma omp barrier


        } //END OF DIJKSTRA WHILE LOOP


/*
        #pragma omp single
        {
            printf("ends\n");
            for(int i=0;i<endsLen;i++){
                printf("%d - %d\n",i,ends[i]);
            }
            printf("settled\n");
            for(int i=0;i<slen;i++){
                printf("%d - %d\n",i,settled[i]);
            }
        }
*/


        /*** 3. DEPENDENCY SUMMATION ***/

        depth = endsLen - 1;
        while(depth > 0){
            start = ends[depth -1];
            end = ends[depth] ;

            #pragma omp for
            for(int i=0; i < end - start; i++){
                int w = settled[start+i];
                //printf("%d being analyzed\n",w);
                int dsw = 0;
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
                    //atomicadd
                    printf("\n**** PARTIAL CONTRIBUTE TO CB, DEP FROM %d PASSING FROM %d = %f  *******\n", root, w, dep[w]);
                    #pragma omp atomic
                        arrayBC[w] = arrayBC[w]+dep[w];
                }
            }
            #pragma omp barrier

            #pragma omp single
            {
                depth = depth - 1;
            }
        } //end of while on DAG


    } //END OF PARALLEL REGION
}


double* BC(int N,int M,int* nodes,Edge* edges /*calNode, calEdge, bcpath.c_str(), ebcpath.c_str(), warp_size*/){
    double* arrayBC;
    arrayBC = malloc(N*sizeof(double));


    printf("I'm Here!\n");

    #pragma omp parallel for
    for(int i=0; i<N; i++){
        //printf("I'm Here! for node %d\n", i);
        processSingleRoot(i,N,M,nodes,edges,arrayBC);
    }

    return arrayBC;
}

