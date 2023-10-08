
#ifndef READ_CDBG_COUNTS_HPP
#define READ_CDBG_COUNTS_HPP



#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Edge.hpp"
#include "Endpoint.hpp"
#include "Build_Params.hpp"
#include "Progress_Tracker.hpp"

#include <cstdint>
#include <string>


template <uint16_t k> class Kmer_SPMC_Iterator;
template <uint16_t k> class Thread_Pool;


// A class to construct compacted read de Bruijn graphs.
template <uint16_t k>
class Read_CdBG_Counts
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (canonical k-mers) of the graph.

    uint64_t edge_count_;    // Number of edges in the underlying graph.
    
    // Members required to keep track of the total number of edges processed across different threads.
    mutable Spin_Lock lock;
    mutable uint64_t edges_processed = 0;
    
    Progress_Tracker progress_tracker;  // Progress tracker for the counts computation task.


    // Distributes the counts computation task — disperses the graph edges (i.e. (k + 1)-mers)
    // parsed by the parser `edge_parser` to the worker threads in the thread pool `thread_pool`,
    // for the edges to be processed by making appropriate state transitions for their endpoints.
    void distribute_counts_computation(Kmer_SPMC_Iterator<k + 1>* edge_parser, Thread_Pool<k>& thread_pool);

    // Processes the edges provided to the thread with id `thread_id` from the parser `edge_parser`,
    // to compute counts.
    void process_edges_counts(Kmer_SPMC_Iterator<k + 1>* edge_parser, uint16_t thread_id);

    // Adds the counts for both k-mers of the edge `e = {u, v}`. Returns `false` iff operation failed
    bool populate_counts(const Edge<k>& e);


public:

    // Consructs a read-CdBG builder object, with the required parameters wrapped in `params`, and uses
    // the Cuttlefish hash table `hash_table`.
    Read_CdBG_Counts(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Computes k-mer counts given the k+1-mer counts in the edge set at path prefix `edge_db_path`.
    void compute_counts(const std::string& edge_db_path);

    
};


template <uint16_t k>
bool Read_CdBG_Counts<k>::populate_counts(const Edge<k>& e)
{
    // Fetch the hash table entry for the vertices associated to the endpoints.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_u = hash_table[e.u().canonical()];
    
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_v = hash_table[e.v().canonical()];

    bucket_u.get_state().set_state(42);
    bucket_v.get_state().set_state(42);


    return hash_table.update_concurrent(bucket_u, bucket_v);
}



#endif
