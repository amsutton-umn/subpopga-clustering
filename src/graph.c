#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <lzma.h>
#include <unistd.h>
#include <error.h>
#include <sys/wait.h>
#include <sys/mman.h>

#include "graph.h"


/*
 * Utility function to check if a file claims to be in xz-compressed
 * format
 */
bool check_if_xz(FILE* file)
{
  int8_t magic[] = {0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00};
  int8_t buf[6];
  int sz = 6;
  sz = fread(buf,1,6,file);
  fseek(file, -sz, SEEK_CUR);
  if (sz == 6){
    if (memcmp(magic,buf,6) == 0) return true;
  }
  return false;
}

void debug(char* s)
{
  fprintf(stderr,"%s\n",s);
}


/*
 * Read a graph from a file and store in data structure
 *
 * The file is in edge list format and can be plain text or
 * xz-compressed
 *
 * Returns the number of edges (if successful), otherwise -1
 *
 */
int read_graph(graph_data* G, FILE* file)
{
  int i,j,k;
  const int edgebuf_size = 1048576;
  int edgebuf[edgebuf_size];
  
  if (check_if_xz(file)){
    debug("Reading compressed file");
    G->m = read_edges_from_compressed(edgebuf,edgebuf_size,file);
  }
  else {
    debug("Reading plaintext file");
    G->m = read_edges_from_plaintext(edgebuf,edgebuf_size,file);
  }
  if (G->m < 1) return -1;

  if (G->m > USHRT_MAX){
    /* Prevent overflow in pair (needed for triangle storage) */
    error(1,ERANGE,"too many edges");
  }

  
  /* determine the number of vertices based on the highest vertex
     number (this could create isolated vertices, but we will just
     ignore these) */
  debug("Determining the vertex count...");
  G->n = 0;
  for (i=0; i<G->m; i++){    
    if (edgebuf[2*i] > G->n) G->n = edgebuf[2*i];
    if (edgebuf[2*i+1] > G->n) G->n = edgebuf[2*i+1];
  }
  G->n++;

  
  /* Allocate memory for graph data structure */
  debug("Allocating memory...");

  G->edge_list    = malloc(G->m*sizeof(pair_t));
  G->adj_list     = malloc(G->n*sizeof(int*));
  G->adj_list_len = malloc(G->n*sizeof(int));

  G->p3_vlist     = malloc(G->n*sizeof(pair_t*));
  G->p3_vlist_len = malloc(G->n*sizeof(int));  
  for (i=0; i<G->n; i++){
    G->adj_list[i] = malloc(G->n*sizeof(int));
    G->p3_vlist[i]  = malloc(G->m*sizeof(pair_t));
    G->adj_list_len[i] = 0;
    G->p3_vlist_len[i]  = 0;
  }

  G->p3_elist     = malloc(G->m*sizeof(int*));
  G->p3_elist_len = malloc(G->m*sizeof(int));
  G->triangle_elist     = malloc(G->m*sizeof(pair_t*));
  G->triangle_elist_len = malloc(G->m*sizeof(int));
  for (i=0; i<G->m; i++){
    G->p3_elist[i]  = malloc(G->m*sizeof(int));
    G->p3_elist_len[i]  = 0;
    G->triangle_elist[i] = malloc(G->m*sizeof(pair_t));
    G->triangle_elist_len[i] = 0;
  }
  
  /* Copy edge buffer to edge list */
  debug("Storing edge list...");
  for (i=0; i<G->m; i++){
    if (edgebuf[2*i] > USHRT_MAX || edgebuf[2*i+1] > USHRT_MAX){
      /* Prevent overflow in pair */
      error(1,ERANGE,"vertex label too large");
    }
    G->edge_list[i] = make_pair(edgebuf[2*i],edgebuf[2*i+1]);
  }

  /* Compute adjacency lists */
  debug("Computing adjacency lists...");
  for (i=0; i<G->m; i++){
    G->adj_list[src(G->edge_list[i])][G->adj_list_len[src(G->edge_list[i])]++] = snk(G->edge_list[i]);
    G->adj_list[snk(G->edge_list[i])][G->adj_list_len[snk(G->edge_list[i])]++] = src(G->edge_list[i]);
    //G->adj_list[edgebuf[2*i]][G->adj_list_len[edgebuf[2*i]]++] = edgebuf[2*i+1];
    //G->adj_list[edgebuf[2*i+1]][G->adj_list_len[edgebuf[2*i+1]]++] = edgebuf[2*i];
  }




  /*  Compute vertex p3 lists */
  debug("Computing and storing vertex P3 lists...");
  for (i=0; i<G->n; i++) {
    for (j=0; j<G->adj_list_len[i]; j++) {
      int v = i;
      int u = G->adj_list[i][j];
      for (k=0; k<G->adj_list_len[u]; k++){
	if (G->adj_list[u][k] != v){
	  int w = G->adj_list[u][k];
	  int triangle=0;
	  int ii;
	  for (ii=0; ii < G->adj_list_len[w]; ii++){
	    if (G->adj_list[w][ii] == v) triangle = 1;
	  }
	  if (!triangle) {
	    G->p3_vlist[v][G->p3_vlist_len[v]++] = make_pair(u,w);
	    if(G->p3_vlist_len[v] > G->m){
	      error(1,ERANGE,"vertex P3 list length");
	    }
	  }
	}
      }	
    }
  }

  /*  Compute edge p3 lists */
  debug("Computing and storing edge P3 and triangle lists...");
  for (i=0; i<G->m; i++) {
    int v = src(G->edge_list[i]);
    int u = snk(G->edge_list[i]);
    int w;
    /* vu is edge i */
    for (w=0; w<G->n; w++){
      int vw = -1;
      int uw = -1;
      if (w == v || w == u) continue;
      for (j=0; j<G->m; j++) {
	if ( pair_eq(G->edge_list[j],v,w) ) vw = j;
	if ( pair_eq(G->edge_list[j],u,w) ) uw = j;
      }      
      /* vw is an edge by uw is not */
      if (vw >= 0 && uw < 0){
	G->p3_elist[i][G->p3_elist_len[i]++] = vw;
      }
      /* uw is an edge but vw is not */
      if (vw < 0 && uw >= 0){
	G->p3_elist[i][G->p3_elist_len[i]++] = uw;
      }
      /* both vw and uw are edges: store triangle */
      if (vw >= 0 && uw >=0){
	G->triangle_elist[i][G->triangle_elist_len[i]++] = make_pair(vw,uw);
      }
    }
  }

  // DEBUG (check this with the picture)
  /* for (i=0; i<G->m; i++) { */
  /*   fprintf(stderr,"Edge %d (%d,%d): ",i,src(G->edge_list[i]),snk(G->edge_list[i])); */
  /*   for (j=0; j<G->p3_elist_len[i]; j++){ */
  /*     fprintf(stderr," (%d,%d) ",src(G->edge_list[G->p3_elist[i][j]]),snk(G->edge_list[G->p3_elist[i][j]])); */
  /*   } */
  /*   fprintf(stderr,"\n"); */
  /* } */

  return G->n;
}

/*
 * Read a list of edges from a plaintext file, storing in edgebuf
 *
 * The bufsize is the size of edgebuf, and cannot be less than twice
 * the number of edges to read.
 *
 * Returns the number of edges read (zero in case of failure)
 */
int read_edges_from_plaintext(int* edgebuf, const int bufsize, FILE* file)
{
  char line[256];
  int ecount = 0;
  int lineno = 0;
  
  while (fgets(line,sizeof(line),file) != NULL) {
    if(line[0] != '#'){
      int u,v,num;
      lineno++;
      num = sscanf(line,"%d %d",&u,&v);
      if (num != 2){
	/* some files contain the number of vertices in the first line */
	fprintf(stderr, "WARNING: line %d contains %d fields.\n",lineno,num);
      }
      else{
	edgebuf[2*ecount] = u;
	edgebuf[2*ecount+1] = v;
	if (2*(++ecount) > bufsize){
	  fprintf(stderr,"ERROR: buffer overflow. Edge buffer size %d is too small.\n",bufsize);
	  return 0;
	}
      }
    }
  }

  return ecount;
}

/*
 * Read edges from a file that has been xz-compressed. The function
 * works by forking a child process. The parent process decompresses
 * the file and writes on a pipe which is read from by the child
 * process by read_edges_from_plaintext.
 *
 * The parent and child must share memory space, which is managed by
 * mmap and munmap.
 */
int read_edges_from_compressed(int* edgebuf, const int bufsize, FILE* file)
{
  /* Map shared memory */
  int* shared_ebuf = mmap(NULL,bufsize*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON,-1,0);
  int* shared_ecount = mmap(NULL,sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON,-1,0);

  int ecount;  
  int fd[2];
  if ( pipe(fd) < 0){
    perror("In read_edges_from_compressed");
    return(0);
  }
  if (fork()==0){
    // read edges plaintext from pipe
    close(fd[1]);
    FILE* readstream = fdopen(fd[0],"r");
    *shared_ecount = read_edges_from_plaintext(shared_ebuf,bufsize,readstream);
    exit(EXIT_SUCCESS);
  }
  else {
    close(fd[0]);
    FILE* writestream = fdopen(fd[1],"w");    
    lzma_stream strm = LZMA_STREAM_INIT;
    lzma_ret ret = lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED);
    if (ret != LZMA_OK) {    
      fprintf(stderr, "Error initializing the lzma decoder\n");
      wait(NULL);
      return 0;
    }
    lzma_action action = LZMA_RUN;
    uint8_t inbuf[BUFSIZ];
    uint8_t outbuf[BUFSIZ];
    bool running = true;
    strm.next_in = NULL;
    strm.avail_in = 0;
    strm.next_out = outbuf;
    strm.avail_out = sizeof(outbuf);

    while (running) {
      if (strm.avail_in == 0 && !feof(file)) {
	strm.next_in = inbuf;
	strm.avail_in = fread(inbuf, 1, sizeof(inbuf),file);
	if (ferror(file)) {
	  fprintf(stderr, "Read error: %s\n",strerror(errno));
	  running = false;
	}
	if (feof(file)) action = LZMA_FINISH;
    }
    ret = lzma_code(&strm, action);

    if (strm.avail_out == 0 || ret == LZMA_STREAM_END) {
      size_t write_size = sizeof(outbuf) - strm.avail_out;
      if (fwrite(outbuf, 1, write_size, writestream)!= write_size) {
	fprintf(stderr, "Write error: %s\n", strerror(errno));
	running = false;
      }	
      strm.next_out = outbuf;
      strm.avail_out = sizeof(outbuf);
    }
      
    if (ret != LZMA_OK) {
      if (ret != LZMA_STREAM_END){
	fprintf(stderr, "Error in lzma decoder: (error code %d)\n", ret);
      }
      running = false;
    }
    }
    fclose(writestream);
    close(fd[1]);    
    lzma_end(&strm);
    wait(NULL);
  }
  /* Clean up */
  ecount = *shared_ecount;
  memcpy(edgebuf,shared_ebuf,sizeof(int)*bufsize);
  munmap(shared_ecount,sizeof(*shared_ecount));
  munmap(shared_ebuf,sizeof(int)*bufsize);
  return ecount;
}
    

void free_graph(graph_data* G)
{
  int i;
  debug("Freeing memory...");
  for (i=0; i<G->n; i++){
    free(G->adj_list[i]);
    free(G->p3_vlist[i]);
  }
  for (i=0; i<G->m; i++){
    free(G->p3_elist[i]);
  }  
  free(G->edge_list);
  free(G->adj_list);
  free(G->adj_list_len);
  free(G->p3_vlist);
  free(G->p3_vlist_len);
  free(G->p3_elist);
  free(G->p3_elist_len);
}









