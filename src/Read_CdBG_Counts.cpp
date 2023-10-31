
#include "Read_CdBG_Counts.hpp"
#include "Edge.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Thread_Pool.hpp"
#include "kseq/kseq.h"
#include "Directed_Kmer.hpp"
#include "fmt/format.h"

#include <chrono>
#include <filesystem>

template <uint16_t k>
Read_CdBG_Counts<k>::Read_CdBG_Counts(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
void Read_CdBG_Counts<k>::compute_counts(const std::string& edge_db_path)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    const Kmer_Container<k + 1> edge_container(edge_db_path);  // Wrapper container for the edge-database.
    Kmer_SPMC_Iterator<k + 1> edge_parser(&edge_container, params.thread_count());  // Parser for the edges from the edge-database.
    edge_count_ = edge_container.size();
    std::cout << "Total number of distinct edges: " << edge_count_ << ".\n";


    const std::string& buckets_file_path = params.buckets_file_path();
    if(!buckets_file_path.empty() && file_exists(buckets_file_path))    // The serialized hash table buckets, saved from some earlier execution, exists.
    {
        std::cout <<    "Found the hash table buckets at file " << buckets_file_path << ".\n"
                        "Loading the buckets.\n";
        hash_table.load_hash_buckets(buckets_file_path);
        std::cout << "Loaded the buckets into memory.\n";
    }
    else
    {
        // Construct a thread pool.
        const uint16_t thread_count = params.thread_count();
        Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::compute_counts_read_space);

        // Launch the reading (and parsing per demand) of the edges from disk.
        edge_parser.launch_production();

        // Launch (multi-threaded) computation of the states.
        const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((edge_count_ / 100.0) / params.thread_count()));
        progress_tracker.setup(edge_count_, thread_load_percentile, "Computing k-mer counts from k+1-mer counts");
        distribute_counts_computation(&edge_parser, thread_pool);

        // Wait for the edges to be depleted from the database.
        edge_parser.seize_production();

        // Wait for the consumer threads to finish parsing and processing the edges.
        thread_pool.close();

        std::cout << "\nNumber of processed edges: " << edges_processed << "\n";


        // Save the hash table buckets, if a file path is provided.
        if(params.save_buckets())
        {
            hash_table.save_hash_buckets(buckets_file_path); 
            std::cout << "Saved the hash buckets at " << buckets_file_path << "\n";
        }
    }

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done computing k-mer counts from KMC edge database. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Counts<k>::distribute_counts_computation(Kmer_SPMC_Iterator<k + 1>* const edge_parser, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_read_dBG_counts_task(edge_parser, idle_thread_id);
    }
}


template <uint16_t k>
void Read_CdBG_Counts<k>::process_edges_counts(Kmer_SPMC_Iterator<k + 1>* const edge_parser, const uint16_t thread_id)
{
    Edge<k> e;  // For the edges to be processed one-by-one; say this is between the vertices `u` and `v`.

    uint64_t edge_count = 0;    // Number of edges processed by this thread.
    uint64_t progress = 0;  // Number of edges processed by the thread; is reset at reaching 1% of its approximate workload.

    while(edge_parser->tasks_expected(thread_id))
        if(edge_parser->value_at(thread_id, e.e()))
        {
            e.configure(hash_table);    // A new edge (k + 1)-mer has been parsed; set information for its two endpoints.

            populate_counts(e);

            edge_count++;
            if(progress_tracker.track_work(++progress))
                progress = 0;
        }
    
    lock.lock();
    edges_processed += edge_count;
    lock.unlock();
}


// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(int, read)

template <uint16_t k>
void Read_CdBG_Counts<k>::write_unitigs_mean_abundances(const std::string& unitigs_path)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    // STEP 2: open the file handler
    FILE* fp = fopen(unitigs_path.c_str(), "r");

    // STEP 3: initialize seq
    kseq_t* parser = kseq_init(fileno(fp));

    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::output_unitigs_with_counts);

    // We're renumbering the unitigs for efficiency, not that the original numbering mattered anyway..
    size_t unitig_id = 0;

    std::string tmp_unitigs_path = unitigs_path+".tmp_abundances";
    output_ = std::ofstream(tmp_unitigs_path);

    // STEP 4: read sequence
    while(kseq_read(parser) >= 0)
    {
        // Multi-threaded
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_output_task(idle_thread_id, parser->seq.s, parser->seq.l, unitig_id, 0);
        thread_pool.wait_completion();
        unitig_id += 1;
    }

    // Close the thread-pool.
    thread_pool.close();

    kseq_destroy(parser);
    fclose(fp);
    output_.close();

    namespace fs = std::filesystem;
    try {
        if (fs::exists(unitigs_path) && fs::is_regular_file(unitigs_path)) {
            fs::remove(unitigs_path);
        }
        fs::rename(tmp_unitigs_path, unitigs_path);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error moving abundance-populated unitigs files: " << e.what() << std::endl;
    }
    
    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Done writing approximate k-mer counts to unitigs. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
void Read_CdBG_Counts<k>::process_unitig_with_counts(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t unitig_id, const size_t dummy)
{
    (void)dummy; (void)thread_id; // not using these values
    const Kmer<k> first_kmer(seq, 0);
    Directed_Kmer<k> kmer(first_kmer);

    float mean_abundance = 0;

    // Scan through the k-mers one-by-one.
    for(size_t kmer_idx = 0; kmer_idx <= seq_len - k; ++kmer_idx)
    {
        Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[kmer.canonical()];
        uint64_t count = decodeAbundance(bucket.get_state().get_state());
        mean_abundance += count;

        if(kmer_idx < seq_len - k)
            kmer.roll_to_next_kmer(seq[kmer_idx + k]);
    }
    mean_abundance /= (seq_len-k+1);
    mean_abundance /= 2; 

    std::string buffer;
    buffer += ">";
    buffer += fmt::format_int(unitig_id).c_str();
    buffer += " ka:f:"; // ka: mean approximate k-mer abundance. 
                        // Highly imprecise due to 2 sources of imprecision: 
                        // - the 8 bit encoding, and also 
                        // - the extraction of k-mer counts from k+1-mer counts
                        // This results in underestimations of low counts
    buffer += fmt::format("{:.1f}", mean_abundance); 
    buffer += "\n";
    buffer += seq;
    buffer += "\n";

    // not the most efficient.. but that'll do
    lock.lock();
    output_ << buffer;
    lock.unlock();
}


// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Counts)
