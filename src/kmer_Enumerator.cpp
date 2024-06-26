
#include "globals.hpp"
#include "kmer_Enumerator.hpp"
#include "Kmer_Container.hpp"
#include "kmer_Enumeration_Stats.hpp"


template <uint16_t k> const std::size_t kmer_Enumerator<k>::min_memory;


template <uint16_t k>
kmer_Enumeration_Stats<k> kmer_Enumerator<k>::enumerate(
    const KMC::InputFileType input_file_type, const std::vector<std::string>& seqs, const uint32_t cutoff, const uint16_t thread_count,
    const std::size_t max_memory, const bool strict_memory, const bool estimate_mem_usage, const double bits_per_kmer,
    const std::string& working_dir_path, const std::string& output_db_path)
{
    // FunnyProgress progress;

    const bool estimate_mem = (k > small_k_threshold && estimate_mem_usage);   // Histogram estimation is not supported for small enough k's yet.
    std::size_t memory = std::max(max_memory, min_memory);
    stage1_params
        .SetInputFileType(input_file_type)
        .SetInputFiles(seqs)
        .SetKmerLen(k)
        .SetNThreads(thread_count)
        .SetTmpPath(working_dir_path)
        .SetEstimateHistogramCfg(estimate_mem ? KMC::EstimateHistogramCfg::ESTIMATE_AND_COUNT_KMERS : KMC::EstimateHistogramCfg::DONT_ESTIMATE)
        // .SetPercentProgressObserver(&progress)
    ;

    if(strict_memory)
        stage1_params
            .SetMaxRamGB(memory)
            .SetSignatureLen(signature_len)
            .SetNBins(bin_count)
        ;

    stage1_results = kmc.RunStage1(stage1_params);


    memory = std::max(
        (estimate_mem ? std::max(memory_limit(solid_kmer_count_approx(cutoff), bits_per_kmer), max_memory) : max_memory),
        min_memory);
    stage2_params
        .SetCutoffMin(cutoff)
        .SetNThreads(thread_count)
        .SetStrictMemoryMode(strict_memory)
//#ifndef CF_VALIDATION_MODE
        .SetCounterMax(counter_max)
//#endif
        .SetOutputFileName(output_db_path)
    ;

    if(strict_memory)
        stage2_params.SetMaxRamGB(memory);

    stage2_results = kmc.RunStage2(stage2_params);
    const std::size_t db_size = Kmer_Container<k>::database_size(output_db_path);

    return kmer_Enumeration_Stats<k>(stage1_results, stage2_results, memory, db_size);
}


template <uint16_t k>
uint64_t kmer_Enumerator<k>::solid_kmer_count_approx(const uint16_t cutoff) const
{
    uint64_t solid_kmer_count = 0;
    for (std::size_t freq = cutoff; freq < stage1_results.estimatedHistogram.size(); ++freq)
        solid_kmer_count += stage1_results.estimatedHistogram[freq];

    return solid_kmer_count;
}


template <uint16_t k>
std::size_t kmer_Enumerator<k>::memory_limit(const uint64_t unique_kmer_count, const double bits_per_kmer) const
{
    const double memory_in_bits = bits_per_kmer * unique_kmer_count;
    const double memory_in_bytes = memory_in_bits / 8.0;
    std::size_t memory_in_gb = static_cast<std::size_t>(memory_in_bytes / (1024 * 1024 * 1024));

    return memory_in_gb;
}



// Template instantiations for the required instances. 
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_ALL, kmer_Enumerator)
