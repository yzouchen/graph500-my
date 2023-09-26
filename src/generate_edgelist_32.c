/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */
/*           Anton Korzh                                                   */


/* These need to be before any possible inclusions of stdint.h or inttypes.h.
 * */

/**
 * This is siplified from main.c to just generate edgelist (V_size=32bit)
*/
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "../generator/make_graph.h"
#include "../generator/utils.h"
#include "aml.h"
#include "common.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

int main(int argc, char** argv) {
	aml_init(&argc,&argv); //includes MPI_Init inside
	setup_globals();

	/* Parse arguments. */
	int SCALE = 16;
	int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
    const char* filename = "Kron.txt";
	if (argc >= 2) SCALE = atoi(argv[1]);
	if (argc >= 3) edgefactor = atoi(argv[2]);
    if (argc >= 4) filename = argv[3];
	if (argc <= 1 || argc >= 5 || SCALE == 0 || edgefactor == 0) {
		if (rank == 0) {
			fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [integer, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	uint64_t seed1 = 2, seed2 = 3;

	const int reuse_file = getenv("REUSEFILE")? 1 : 0;
	/* If filename is NULL, store data in memory */

	tuple_graph tg;
	tg.nglobaledges = (int64_t)(edgefactor) << SCALE;
	int64_t nglobalverts = (int64_t)(1) << SCALE;

	tg.data_in_file = 0;
	tg.write_file = 1;

	/* Make the raw graph edges. */
	double make_graph_start = MPI_Wtime();
	if( !tg.data_in_file || tg.write_file )
	{
		/* Spread the two 64-bit numbers into five nonzero values in the correct
		 * range. */
		uint_fast32_t seed[5];
		make_mrg_seed(seed1, seed2, seed);

		/* As the graph is being generated, also keep a bitmap of vertices with
		 * incident edges.  We keep a grid of processes, each row of which has a
		 * separate copy of the bitmap (distributed among the processes in the
		 * row), and then do an allreduce at the end.  This scheme is used to avoid
		 * non-local communication and reading the file separately just to find BFS
		 * roots. */
		MPI_Offset nchunks_in_file = (tg.nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE;
		int64_t bitmap_size_in_bytes = int64_min(BITMAPSIZE, (nglobalverts + CHAR_BIT - 1) / CHAR_BIT);
		if (bitmap_size_in_bytes * size * CHAR_BIT < nglobalverts) {
			bitmap_size_in_bytes = (nglobalverts + size * CHAR_BIT - 1) / (size * CHAR_BIT);
		}
		int ranks_per_row = tg.data_in_file ? ((nglobalverts + CHAR_BIT - 1) / CHAR_BIT + bitmap_size_in_bytes - 1) / bitmap_size_in_bytes : 1;
		int nrows = size / ranks_per_row;
		int my_row = -1, my_col = -1;
		MPI_Comm cart_comm;
		{
			int dims[2] = {size / ranks_per_row, ranks_per_row};
			int periods[2] = {0, 0};
			MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
		}
		int in_generating_rectangle = 0;
		if (cart_comm != MPI_COMM_NULL) {
			in_generating_rectangle = 1;
			{
				int dims[2], periods[2], coords[2];
				MPI_Cart_get(cart_comm, 2, dims, periods, coords);
				my_row = coords[0];
				my_col = coords[1];
			}
			MPI_Comm this_col;
			MPI_Comm_split(cart_comm, my_col, my_row, &this_col);
			MPI_Comm_free(&cart_comm);
			/* Every rank in a given row creates the same vertices (for updating the
			 * bitmap); only one writes them to the file (or final memory buffer). */
			packed_edge* buf = (packed_edge*)xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge));

			MPI_Offset block_limit = (nchunks_in_file + nrows - 1) / nrows;
			/* fprintf(stderr, "%d: nchunks_in_file = %" PRId64 ", block_limit = %" PRId64 " in grid of %d rows, %d cols\n", rank, (int64_t)nchunks_in_file, (int64_t)block_limit, nrows, ranks_per_row); */
			if (tg.data_in_file) {
				tg.edgememory_size = 0;
				tg.edgememory = NULL;
			} else {
				int my_pos = my_row + my_col * nrows;
				int last_pos = (tg.nglobaledges % ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row) != 0) ?
					(tg.nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row) :
					-1;
				int64_t edges_left = tg.nglobaledges % FILE_CHUNKSIZE;
				int64_t nedges = FILE_CHUNKSIZE * (tg.nglobaledges / ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row)) +
					FILE_CHUNKSIZE * (my_pos < (tg.nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row)) +
					(my_pos == last_pos ? edges_left : 0);
				/* fprintf(stderr, "%d: nedges = %" PRId64 " of %" PRId64 "\n", rank, (int64_t)nedges, (int64_t)tg.nglobaledges); */
				tg.edgememory_size = nedges;
			}

            /* ADD, write file */
            FILE *fp;
            fp = fopen(filename, "w");

			MPI_Offset block_idx;
			for (block_idx = 0; block_idx < block_limit; ++block_idx) {
				/* fprintf(stderr, "%d: On block %d of %d\n", rank, (int)block_idx, (int)block_limit); */
				MPI_Offset start_edge_index = int64_min(FILE_CHUNKSIZE * (block_idx * nrows + my_row), tg.nglobaledges);
				MPI_Offset edge_count = int64_min(tg.nglobaledges - start_edge_index, FILE_CHUNKSIZE);
				// packed_edge* actual_buf = (!tg.data_in_file && block_idx % ranks_per_row == my_col) ?
				// 	tg.edgememory + FILE_CHUNKSIZE * (block_idx / ranks_per_row) :
				// 	buf;

                packed_edge* actual_buf = buf;

				/* fprintf(stderr, "%d: My range is [%" PRId64 ", %" PRId64 ") %swriting into index %" PRId64 "\n", rank, (int64_t)start_edge_index, (int64_t)(start_edge_index + edge_count), (my_col == (block_idx % ranks_per_row)) ? "" : "not ", (int64_t)(FILE_CHUNKSIZE * (block_idx / ranks_per_row))); */
				// if (!tg.data_in_file && block_idx % ranks_per_row == my_col) {
				// 	assert (FILE_CHUNKSIZE * (block_idx / ranks_per_row) + edge_count <= tg.edgememory_size);
				// }
				if (tg.write_file) {
					generate_kronecker_range(seed, SCALE, start_edge_index, start_edge_index + edge_count, actual_buf);
					
                    int i = 0;
                    for(i = 0; i < edge_count; i++){
                        fprintf(fp, "%d %d\n", actual_buf[i].v0, actual_buf[i].v1);
                    }
                    
                    if (tg.data_in_file && my_col == (block_idx % ranks_per_row)) { /* Try to spread writes among ranks */
						MPI_File_write_at(tg.edgefile, start_edge_index, actual_buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
					}
				} else {
					/* All read rather than syncing up for a row broadcast. */
					MPI_File_read_at(tg.edgefile, start_edge_index, actual_buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
				}
			}
            fclose(fp);
			free(buf);
			MPI_Comm_free(&this_col);
		} else {
			tg.edgememory = NULL;
			tg.edgememory_size = 0;
		}
		// MPI_Allreduce(&tg.edgememory_size, &tg.max_edgememory_size, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
		if (tg.data_in_file && tg.write_file) {
			MPI_File_sync(tg.edgefile);
		}
	}

	double make_graph_stop = MPI_Wtime();
	double make_graph_time = make_graph_stop - make_graph_start;
	if (rank == 0) { /* Not an official part of the results */
		fprintf(stderr, "graph_generation:               %f s\n", make_graph_time);
	}

	/* Make user's graph data structure. */
	// double data_struct_start = MPI_Wtime();
	// make_graph_data_structure(&tg);
	// double data_struct_stop = MPI_Wtime();
	// double data_struct_time = data_struct_stop - data_struct_start;
	// if (rank == 0) { /* Not an official part of the results */
	// 	fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
	// }

	return 0;
}
