
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
#include <fstream>

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

    std::ofstream output_;

public:

    // Consructs a read-CdBG builder object, with the required parameters wrapped in `params`, and uses
    // the Cuttlefish hash table `hash_table`.
    Read_CdBG_Counts(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Computes k-mer counts given the k+1-mer counts in the edge set at path prefix `edge_db_path`.
    void compute_counts(const std::string& edge_db_path);

    // Write abundances to the unitigs FASTA file.
    void write_unitigs_mean_abundances(const std::string& unitigs_path);

    // Process all unitigs and compute mean average count, writing the result to disk, replacing the original FASTA file
    void process_unitig_with_counts(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t dummy1, const size_t dummy2);

    // Encode and decode abundances stored in 8 bits using a lossy encoding
    uint8_t encodeAbundance(uint32_t abundance);
    uint32_t decodeAbundance(uint8_t encoded);

};

// inspired by Guillaum's trick in GATB
// https://github.com/GATB/gatb-core/blob/0a9ffd598259eb13a1257ffb28236f1510dfbefc/gatb-core/src/gatb/tools/collections/impl/MapMPHF.hpp#L84
// Convert a normal abundance value to its 8-bit encoded form
template <uint16_t k>
uint8_t Read_CdBG_Counts<k>::encodeAbundance(uint32_t abundance) {
    if (abundance <= 70) {
        return static_cast<uint8_t>(abundance);
    } else if (abundance <= 100) {
        return 70 + (abundance - 70) / 2;
    } else if (abundance <= 500) {
        return 85 + (abundance - 100) / 10;
    } else if (abundance <= 1000) {
        return 125 + (abundance - 500) / 20;
    } else if (abundance <= 5000) {
        return 150 + (abundance - 1000) / 100;
    } else if (abundance <= 10000) {
        return 190 + (abundance - 5000) / 200;
    } else if (abundance <= 50000) {
        return 215 + (abundance - 10000) / 1000;
    } else {
        // Handle values above 50000 if necessary
        return 255;
    }
}

// Decode an 8-bit encoded abundance value back to its original form
template <uint16_t k>
uint32_t Read_CdBG_Counts<k>::decodeAbundance(uint8_t encoded) {
    if (encoded <= 70) {
        return static_cast<uint32_t>(encoded);
    } else if (encoded <= 85) {
        return 70 + (encoded - 70) * 2;
    } else if (encoded <= 125) {
        return 100 + (encoded - 85) * 10;
    } else if (encoded <= 150) {
        return 500 + (encoded - 125) * 20;
    } else if (encoded <= 190) {
        return 1000 + (encoded - 150) * 100;
    } else if (encoded <= 215) {
        return 5000 + (encoded - 190) * 200;
    } else {
        return 10000 + (encoded - 215) * 1000;
    }
}

template <uint16_t k>
bool Read_CdBG_Counts<k>::populate_counts(const Edge<k>& e)
{
    // Fetch the hash table entry for the vertices associated to the endpoints.

    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_u = hash_table[e.u().canonical()];
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket_v = hash_table[e.v().canonical()];

    uint32_t abundance_u = (uint32_t)(e.abundance()) + decodeAbundance(bucket_u.get_state().get_state());
    uint32_t abundance_v = (uint32_t)(e.abundance()) + decodeAbundance(bucket_v.get_state().get_state());
    
    bucket_u.get_state().set_state(encodeAbundance(abundance_u));
    bucket_v.get_state().set_state(encodeAbundance(abundance_v));

    return hash_table.update_concurrent(bucket_u, bucket_v);
}



#endif
