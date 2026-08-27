#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "../ramses.h"
#include "fof.h"
#include "distributed.h"

#define LABEL_RANK_SHIFT 48
#define LABEL_LOCAL_MASK ((UINT64_C(1) << LABEL_RANK_SHIFT) - UINT64_C(1))
#define MAX_LABEL_ROUNDS_FACTOR 4

typedef struct FaceParticle {
	particle p;
	size_t component;
} FaceParticle;

typedef struct ComponentEdge {
	size_t local;
	size_t remote;
} ComponentEdge;

typedef struct BoundaryNode {
	size_t group;
	size_t component;
	int side;
} BoundaryNode;

typedef struct LabelUpdate {
	size_t target;
	uint64_t label;
} LabelUpdate;

typedef struct MemberRecord {
	uint64_t label;
	particle p;
} MemberRecord;

static void *checked_malloc(size_t bytes)
{
	void *value = malloc(bytes > 0 ? bytes : 1);
	if(value == NULL){
		fprintf(stderr,"OPFoF distributed merge: cannot allocate %zu bytes\n",bytes);
		MPI_Abort(MPI_COMM_WORLD,92);
	}
	return value;
}

static void *checked_realloc(void *old, size_t bytes)
{
	void *value = realloc(old,bytes > 0 ? bytes : 1);
	if(value == NULL){
		fprintf(stderr,"OPFoF distributed merge: cannot reallocate %zu bytes\n",bytes);
		MPI_Abort(MPI_COMM_WORLD,92);
	}
	return value;
}

static void big_sendrecv_bytes(
		const void *send, size_t send_bytes, int dest, int send_tag,
		void *recv, size_t recv_bytes, int src, int recv_tag)
{
	const size_t chunk_limit = 500000000;
	size_t send_chunks = (send_bytes+chunk_limit-1)/chunk_limit;
	size_t recv_chunks = (recv_bytes+chunk_limit-1)/chunk_limit;
	size_t nrequest = send_chunks+recv_chunks, at = 0, offset;
	MPI_Request *request = checked_malloc(sizeof(*request)*nrequest);
	for(offset=0;offset<recv_bytes;offset+=chunk_limit){
		size_t amount = recv_bytes-offset;
		if(amount > chunk_limit) amount = chunk_limit;
		MPI_Irecv((char*)recv+offset,(int)amount,MPI_BYTE,src,recv_tag,
				MPI_COMM_WORLD,&request[at++]);
	}
	for(offset=0;offset<send_bytes;offset+=chunk_limit){
		size_t amount = send_bytes-offset;
		if(amount > chunk_limit) amount = chunk_limit;
		MPI_Isend((char*)send+offset,(int)amount,MPI_BYTE,dest,send_tag,
				MPI_COMM_WORLD,&request[at++]);
	}
	if(nrequest > 0) MPI_Waitall((int)nrequest,request,MPI_STATUSES_IGNORE);
	free(request);
}

static int compare_edge(const void *aa, const void *bb)
{
	const ComponentEdge *a = aa, *b = bb;
	if(a->local < b->local) return -1;
	if(a->local > b->local) return 1;
	if(a->remote < b->remote) return -1;
	if(a->remote > b->remote) return 1;
	return 0;
}

static int compare_boundary_node(const void *aa, const void *bb)
{
	const BoundaryNode *a = aa, *b = bb;
	if(a->group < b->group) return -1;
	if(a->group > b->group) return 1;
	if(a->side < b->side) return -1;
	if(a->side > b->side) return 1;
	if(a->component < b->component) return -1;
	if(a->component > b->component) return 1;
	return 0;
}

static int compare_member(const void *aa, const void *bb)
{
	const MemberRecord *a = aa, *b = bb;
	if(a->label < b->label) return -1;
	if(a->label > b->label) return 1;
	return 0;
}

static MemberRecord *exchange_members_left(
		const MemberRecord *send, size_t nsend, size_t *nrecv,
		int left, int right, int hop)
{
	const size_t byte_limit = 500000000;
	size_t send_bytes = nsend*sizeof(*send), recv_bytes, offset;
	size_t send_chunks, recv_chunks, nrequest, request_index = 0;
	MemberRecord *recv;
	MPI_Request *request;
	MPI_Status status;
	MPI_Sendrecv(&nsend,sizeof(size_t),MPI_BYTE,left,501+hop,
			nrecv,sizeof(size_t),MPI_BYTE,right,501+hop,MPI_COMM_WORLD,&status);
	recv = checked_malloc(sizeof(*recv)*(*nrecv));
	recv_bytes = (*nrecv)*sizeof(*recv);
	send_chunks = (send_bytes+byte_limit-1)/byte_limit;
	recv_chunks = (recv_bytes+byte_limit-1)/byte_limit;
	nrequest = send_chunks+recv_chunks;
	request = checked_malloc(sizeof(*request)*nrequest);
	/* Post every receive before any send.  Zero-byte peers need no payload
	 * request because the particle counts have already matched above. */
	for(offset=0;offset<recv_bytes;offset+=byte_limit){
		size_t amount = recv_bytes-offset;
		if(amount > byte_limit) amount = byte_limit;
		MPI_Irecv((char*)recv+offset,(int)amount,MPI_BYTE,right,601+hop,
				MPI_COMM_WORLD,&request[request_index++]);
	}
	for(offset=0;offset<send_bytes;offset+=byte_limit){
		size_t amount = send_bytes-offset;
		if(amount > byte_limit) amount = byte_limit;
		MPI_Isend((char*)send+offset,(int)amount,MPI_BYTE,left,601+hop,
				MPI_COMM_WORLD,&request[request_index++]);
	}
	if(nrequest > 0) MPI_Waitall((int)nrequest,request,MPI_STATUSES_IGNORE);
	free(request);
	return recv;
}

static size_t unique_edges(ComponentEdge *edge, size_t count)
{
	size_t i, out;
	if(count == 0) return 0;
	qsort(edge,count,sizeof(*edge),compare_edge);
	out = 1;
	for(i=1;i<count;i++){
		if(compare_edge(&edge[i],&edge[out-1]) != 0) edge[out++] = edge[i];
	}
	return out;
}

static void ensure_tree(FoFTStruct **tree, long long *capacity, size_t particles)
{
	long long needed = (long long)(particles/NODE_HAVE_PARTICLE*10 + 1024);
	if(needed > *capacity){
		*capacity = needed;
		*tree = checked_realloc(*tree,sizeof(**tree)*(size_t)(*capacity));
	}
}

static ComponentEdge *find_right_edges(
		FoFTPtlStruct *ptl, size_t np, FoFTStruct **tree,
		long long *tree_capacity, Box box, POSTYPE boundary_link,
		POSTYPE domain_low, POSTYPE domain_high, POSTYPE input_fof_link,
		int periodic, POSTYPE box_x, POSTYPE box_y, POSTYPE box_z,
		size_t *edge_count)
{
	int rank, nrank, left, right;
	size_t i, nlower = 0, nupper = 0, nremote = 0, ntotal;
	FaceParticle *lower, *upper, *remote;
	FoFTPtlStruct *work, *linked;
	size_t *component;
	int *side;
	BoundaryNode *nodes;
	ComponentEdge *edges = NULL;
	size_t nedge = 0, edge_capacity = 0, ngroup = 0;
	MPI_Status status;
	particle seed;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	left = (rank-1+nrank)%nrank;
	right = (rank+1)%nrank;
	for(i=0;i<np;i++){
		if(ptl[i].z <= domain_low+boundary_link) nlower++;
		if(ptl[i].z >= domain_high-boundary_link) nupper++;
	}
	lower = checked_malloc(sizeof(*lower)*nlower);
	upper = checked_malloc(sizeof(*upper)*nupper);
	nlower = nupper = 0;
	for(i=0;i<np;i++){
		if(ptl[i].z <= domain_low+boundary_link){
			lower[nlower].p = ptl[i];
			lower[nlower++].component = ptl[i].haloindx;
		}
		if(ptl[i].z >= domain_high-boundary_link){
			upper[nupper].p = ptl[i];
			upper[nupper++].component = ptl[i].haloindx;
		}
	}

	MPI_Sendrecv(&nlower,sizeof(size_t),MPI_BYTE,left,401,
			&nremote,sizeof(size_t),MPI_BYTE,right,401,MPI_COMM_WORLD,&status);
	remote = checked_malloc(sizeof(*remote)*nremote);
	big_sendrecv_bytes(lower,nlower*sizeof(*lower),left,402,
			remote,nremote*sizeof(*remote),right,402);
	free(lower);

	ntotal = nupper+nremote;
	if(ntotal == 0){
		free(upper);
		free(remote);
		*edge_count = 0;
		return NULL;
	}
	work = checked_malloc(sizeof(*work)*ntotal);
	linked = checked_malloc(sizeof(*linked)*ntotal);
	component = checked_malloc(sizeof(*component)*ntotal);
	side = checked_malloc(sizeof(*side)*ntotal);
	for(i=0;i<nupper;i++){
		work[i] = upper[i].p;
		component[i] = upper[i].component;
		side[i] = 0;
	}
	for(i=0;i<nremote;i++){
		work[nupper+i] = remote[i].p;
		component[nupper+i] = remote[i].component;
		side[nupper+i] = 1;
	}
	free(upper);
	free(remote);
	for(i=0;i<ntotal;i++){
		work[i].sibling = (i+1 < ntotal) ? &work[i+1] : NULL;
		work[i].haloindx = (size_t)-1;
	}
	ensure_tree(tree,tree_capacity,ntotal);
	FoF_Make_Tree(*tree,work,ntotal,box);
	for(i=0;i<ntotal;i++){
		if(work[i].included == NO){
			seed.x = work[i].x;
			seed.y = work[i].y;
			seed.z = work[i].z;
			seed.link02 = work[i].link02;
			if(periodic) pnew_fof_link(&seed,input_fof_link,*tree,work,
					linked,ngroup,box_x,box_y,box_z);
			else new_fof_link(&seed,input_fof_link,*tree,work,linked,ngroup);
			ngroup++;
		}
	}
	free(linked);
	nodes = checked_malloc(sizeof(*nodes)*ntotal);
	for(i=0;i<ntotal;i++){
		nodes[i].group = work[i].haloindx;
		nodes[i].component = component[i];
		nodes[i].side = side[i];
	}
	free(work);
	free(component);
	free(side);
	qsort(nodes,ntotal,sizeof(*nodes),compare_boundary_node);

	for(i=0;i<ntotal;){
		size_t end = i+1, local_begin, local_end, remote_begin, remote_end, a, b;
		while(end < ntotal && nodes[end].group == nodes[i].group) end++;
		local_begin = i;
		while(local_begin < end && nodes[local_begin].side != 0) local_begin++;
		local_end = local_begin;
		while(local_end < end && nodes[local_end].side == 0) local_end++;
		remote_begin = local_end;
		while(remote_begin < end && nodes[remote_begin].side != 1) remote_begin++;
		remote_end = remote_begin;
		while(remote_end < end && nodes[remote_end].side == 1) remote_end++;
		for(a=local_begin;a<local_end;a++){
			if(a > local_begin && nodes[a].component == nodes[a-1].component) continue;
			for(b=remote_begin;b<remote_end;b++){
				if(b > remote_begin && nodes[b].component == nodes[b-1].component) continue;
				if(nedge == edge_capacity){
					edge_capacity = edge_capacity ? 2*edge_capacity : 256;
					edges = checked_realloc(edges,sizeof(*edges)*edge_capacity);
				}
				edges[nedge].local = nodes[a].component;
				edges[nedge++].remote = nodes[b].component;
			}
		}
		i = end;
	}
	free(nodes);
	nedge = unique_edges(edges,nedge);
	*edge_count = nedge;
	printf("P%d found %zu component links across its right interface\n",rank,nedge);
	return edges;
}

static ComponentEdge *mirror_left_edges(
		ComponentEdge *right_edges, size_t nright, size_t *nleft)
{
	int rank, nrank, left, right;
	size_t i;
	ComponentEdge *incoming, *left_edges;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	left = (rank-1+nrank)%nrank;
	right = (rank+1)%nrank;
	MPI_Sendrecv(&nright,sizeof(size_t),MPI_BYTE,right,411,
			nleft,sizeof(size_t),MPI_BYTE,left,411,MPI_COMM_WORLD,&status);
	incoming = checked_malloc(sizeof(*incoming)*(*nleft));
	big_sendrecv_bytes(right_edges,nright*sizeof(*right_edges),right,412,
			incoming,(*nleft)*sizeof(*incoming),left,412);
	left_edges = checked_malloc(sizeof(*left_edges)*(*nleft));
	for(i=0;i<*nleft;i++){
		left_edges[i].local = incoming[i].remote;
		left_edges[i].remote = incoming[i].local;
	}
	free(incoming);
	*nleft = unique_edges(left_edges,*nleft);
	return left_edges;
}

static uint64_t *converge_component_labels(
		size_t nhalo, ComponentEdge *right_edges, size_t nright,
		ComponentEdge *left_edges, size_t nleft)
{
	int rank, nrank, left, right, round;
	size_t i;
	uint64_t *label;
	LabelUpdate *send_right, *send_left, *recv_left, *recv_right;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	left = (rank-1+nrank)%nrank;
	right = (rank+1)%nrank;
	if((uint64_t)rank >= (UINT64_C(1) << (64-LABEL_RANK_SHIFT)) ||
			(uint64_t)nhalo > LABEL_LOCAL_MASK){
		fprintf(stderr,"P%d component label range overflow\n",rank);
		MPI_Abort(MPI_COMM_WORLD,93);
	}
	label = checked_malloc(sizeof(*label)*nhalo);
	for(i=0;i<nhalo;i++) label[i] = ((uint64_t)rank << LABEL_RANK_SHIFT) | i;
	send_right = checked_malloc(sizeof(*send_right)*nright);
	send_left = checked_malloc(sizeof(*send_left)*nleft);
	recv_left = checked_malloc(sizeof(*recv_left)*nleft);
	recv_right = checked_malloc(sizeof(*recv_right)*nright);
	for(round=0;round<MAX_LABEL_ROUNDS_FACTOR*nrank+16;round++){
		int changed = 0, global_changed = 0;
		for(i=0;i<nright;i++){
			send_right[i].target = right_edges[i].remote;
			send_right[i].label = label[right_edges[i].local];
		}
		for(i=0;i<nleft;i++){
			send_left[i].target = left_edges[i].remote;
			send_left[i].label = label[left_edges[i].local];
		}
		big_sendrecv_bytes(send_right,nright*sizeof(*send_right),right,421,
				recv_left,nleft*sizeof(*recv_left),left,421);
		big_sendrecv_bytes(send_left,nleft*sizeof(*send_left),left,422,
				recv_right,nright*sizeof(*recv_right),right,422);
		for(i=0;i<nleft;i++){
			if(recv_left[i].target >= nhalo) MPI_Abort(MPI_COMM_WORLD,94);
			if(recv_left[i].label < label[recv_left[i].target]){
				label[recv_left[i].target] = recv_left[i].label;
				changed = 1;
			}
		}
		for(i=0;i<nright;i++){
			if(recv_right[i].target >= nhalo) MPI_Abort(MPI_COMM_WORLD,94);
			if(recv_right[i].label < label[recv_right[i].target]){
				label[recv_right[i].target] = recv_right[i].label;
				changed = 1;
			}
		}
		MPI_Allreduce(&changed,&global_changed,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
		if(!global_changed){
			if(rank==0) printf("Distributed component labels converged in %d rounds\n",round+1);
			break;
		}
	}
	if(round == MAX_LABEL_ROUNDS_FACTOR*nrank+16){
		fprintf(stderr,"P%d component labels did not converge\n",rank);
		MPI_Abort(MPI_COMM_WORLD,95);
	}
	free(send_right);
	free(send_left);
	free(recv_left);
	free(recv_right);
	return label;
}

static void route_and_write(
		FoFTPtlStruct *ptl, size_t np, uint64_t *label,
		const unsigned char *keep_component,
		const char *halofile, const char *memberfile)
{
	int rank, nrank, hop, peer, left, right;
	size_t i, selected = 0, owned_count = 0, owned_capacity = 1;
	size_t transit_count = 0;
	MemberRecord *owned, *transit;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	left = (rank-1+nrank)%nrank;
	right = (rank+1)%nrank;
	for(i=0;i<np;i++){
		if(!keep_component[ptl[i].haloindx]) continue;
		int owner = (int)(label[ptl[i].haloindx] >> LABEL_RANK_SHIFT);
		if(owner < 0 || owner >= nrank) MPI_Abort(MPI_COMM_WORLD,97);
		selected++;
		if(owner == rank) owned_count++;
		else transit_count++;
	}
	owned_capacity = owned_count > 0 ? owned_count : 1;
	owned = checked_malloc(sizeof(*owned)*owned_capacity);
	transit = checked_malloc(sizeof(*transit)*transit_count);
	{
		size_t owned_at = 0, transit_at = 0;
		for(i=0;i<np;i++){
			uint64_t value;
			int owner;
			MemberRecord record;
			if(!keep_component[ptl[i].haloindx]) continue;
			value = label[ptl[i].haloindx];
			owner = (int)(value >> LABEL_RANK_SHIFT);
			record.label = value;
			record.p = ptl[i];
			if(owner == rank) owned[owned_at++] = record;
			else transit[transit_at++] = record;
		}
	}
	printf("P%d retains %zu of %zu particles; %zu start in the neighbour ring\n",
			rank,selected,np,transit_count); fflush(stdout);
	for(hop=0;hop<nrank-1;hop++){
		size_t incoming_count = 0, next_count = 0;
		MemberRecord *incoming = exchange_members_left(
				transit,transit_count,&incoming_count,left,right,hop);
		MemberRecord *next = checked_malloc(sizeof(*next)*incoming_count);
		free(transit);
		for(i=0;i<incoming_count;i++){
			int owner = (int)(incoming[i].label >> LABEL_RANK_SHIFT);
			if(owner == rank){
				if(owned_count == owned_capacity){
					owned_capacity = 2*owned_capacity+1024;
					owned = checked_realloc(owned,sizeof(*owned)*owned_capacity);
				}
				owned[owned_count++] = incoming[i];
			}
			else next[next_count++] = incoming[i];
		}
		free(incoming);
		transit = next;
		transit_count = next_count;
	}
	{
		unsigned long long local_left = transit_count, global_left = 0;
		MPI_Allreduce(&local_left,&global_left,1,MPI_UNSIGNED_LONG_LONG,
				MPI_SUM,MPI_COMM_WORLD);
		if(global_left != 0){
			if(rank==0) fprintf(stderr,"%llu particles did not reach their ring owner\n",
					global_left);
			MPI_Abort(MPI_COMM_WORLD,98);
		}
	}
	free(transit);
	printf("P%d completed adjacent-ring owner routing (%zu owned particles)\n",
			rank,owned_count); fflush(stdout);
	printf("P%d starts sorting %zu owned particles\n",rank,owned_count); fflush(stdout);
	qsort(owned,owned_count,sizeof(*owned),compare_member);
	printf("P%d completed owned-particle sorting\n",rank); fflush(stdout);
	{
		size_t owned_halo = 0;
		FoFTPtlStruct *owned_particle = checked_malloc(sizeof(*owned_particle)*owned_count);
		HaloBound *bound;
		uint64_t previous = UINT64_MAX;
		for(i=0;i<owned_count;i++){
			if(i==0 || owned[i].label != previous){
				previous = owned[i].label;
				owned_halo++;
			}
			owned_particle[i] = owned[i].p;
			owned_particle[i].haloindx = owned_halo-1;
		}
		free(owned);
		bound = checked_malloc(sizeof(*bound)*owned_halo);
		for(peer=0;peer<nrank;peer++){
			if(peer==rank && owned_count > 0){
				WriteAllHalo(owned_halo,bound,owned_particle,owned_count,NULL,
						(char*)halofile,(char*)memberfile);
				printf("P%d wrote %zu globally owned components (%zu particles)\n",
						rank,owned_halo,owned_count);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		free(bound);
		free(owned_particle);
	}
}

void DistributedMergeAndWrite(
		FoFTPtlStruct *ptl, size_t np, size_t nhalo,
		FoFTStruct **tree, long long *tree_capacity, Box box,
		POSTYPE boundary_link, POSTYPE domain_low, POSTYPE domain_high,
		POSTYPE input_fof_link, int periodic,
		POSTYPE box_x, POSTYPE box_y, POSTYPE box_z,
		const char *halofile, const char *memberfile)
{
	size_t nright = 0, nleft = 0;
	size_t i;
	ComponentEdge *right_edges = find_right_edges(
			ptl,np,tree,tree_capacity,box,boundary_link,domain_low,domain_high,
			input_fof_link,periodic,box_x,box_y,box_z,&nright);
	ComponentEdge *left_edges = mirror_left_edges(right_edges,nright,&nleft);
	uint64_t *label = converge_component_labels(
			nhalo,right_edges,nright,left_edges,nleft);
	unsigned char *keep = calloc(nhalo,sizeof(*keep));
	size_t *local_size = calloc(nhalo,sizeof(*local_size));
	if(keep == NULL || local_size == NULL) MPI_Abort(MPI_COMM_WORLD,96);
	for(i=0;i<np;i++) local_size[ptl[i].haloindx]++;
	for(i=0;i<nhalo;i++){
		if(local_size[i] >= MinNumMem) keep[i] = 1;
	}
	/* Every endpoint of the distributed component graph must be retained;
	 * several individually unresolved pieces may form a resolved halo. */
	for(i=0;i<nright;i++) keep[right_edges[i].local] = 1;
	for(i=0;i<nleft;i++) keep[left_edges[i].local] = 1;
	free(right_edges);
	free(left_edges);
	free(local_size);
	route_and_write(ptl,np,label,keep,halofile,memberfile);
	free(keep);
	free(label);
}
