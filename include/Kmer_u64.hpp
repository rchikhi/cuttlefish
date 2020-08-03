
#ifndef KMER_U64_HPP
#define KMER_U64_HPP



#include "Kmer_Utility.hpp"
#include "kmc_api/kmc_file.h"

#include <string>
#include <iostream>


class Kmer_Hasher;


class Kmer_u64: public Kmer_Utility
{
    friend class Kmer_Hasher;

private:

    static uint16_t k;  // The k-parameter.
    uint64_t kmer = 0;  // The 64-bit encoding of the underlying k-mer.
    static uint64_t bitmask_MSN;    // Bitmask used to clear the most significant nucleotide character, i.e. the first nucleotide of the k-mer which is at the bits `2k-1 : 2k-2`.


    // Returns the 64-bit encoding of the k-mer.
    uint64_t to_u64() const;
    

public:

    // Default constructs the k-mer with a 0-value, equivalent to "AA...A".
    Kmer_u64(): kmer(0)
    {}

    // Constructs a k-mer from the provided string `label`.
    Kmer_u64(const std::string& label);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer_u64(const char* label, const size_t kmer_idx);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer_u64(const std::string& label, const size_t kmer_idx);

    // Constructs a k-mer from `kmer_api` which is a k-mer object built from KMC.
    Kmer_u64(const CKmerAPI& kmer_api);

    // Copy constructs the k-mer from another k-mer `rhs`.
    Kmer_u64(const Kmer_u64& rhs);

    // Copy assignment operator.
    Kmer_u64& operator=(const Kmer_u64& rhs) = default;

    // Sets the value of the `k` parameter across the `Kmer_u64` class.
    static void set_k(const uint16_t k);

    // Returns the reverese complement of the k-mer.
    Kmer_u64 reverse_complement() const;

    // Returns true iff the bitwise encoding of this k-mer is lesser to the
    // encoding of the other k-mer `rhs`.
    bool operator<(const Kmer_u64& rhs) const;

    // Returns true iff this k-mer is identical to the other k-mer `rhs`.
    bool operator==(const Kmer_u64& rhs) const;

    // Returns `true` iff the k-mer is in the forward direction relative to
    // the other k-mer `kmer_hat`.
    bool in_forward(const Kmer_u64& kmer_hat) const;

    // Transforms this k-mer by chopping off the first nucleotide and
    // appending the next nucleotide `next_nucl` to the end, i.e.
    // rolls the k-mer by one nucleotide. Also sets the passed reverse
    // complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(const char next_nucl, Kmer_u64& rev_compl);

    // Returns the canonical version of the k-mer, comparing it to its
    // reverse complement `rev_compl`.
    Kmer_u64 canonical(const Kmer_u64& rev_compl) const;

    // Returns the canonical version of the k-mer.
    Kmer_u64 canonical() const;

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Gets the k-mer from the KMC api object `kmer_api`.
    void from_CKmerAPI(const CKmerAPI& kmer_api);

    // Returns the k-parameter.
    static uint16_t get_k();
    
    // For debugging purposes.
    friend std::ostream& operator<<(std::ostream& out, const Kmer_u64& kmer);
};


inline Kmer_u64::Kmer_u64(const char* label, const size_t kmer_idx)
{
    kmer = 0;

    for(size_t idx = kmer_idx; idx < kmer_idx + k; ++idx)
    {
        uint8_t nucleotide = map_nucleotide(label[idx]);
        kmer = (kmer << 2) | nucleotide;
    }
}


inline Kmer_u64::Kmer_u64(const CKmerAPI& kmer_api)
{
    kmer = 0;

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        // uint8_t nucleotide = map_nucleotide(kmer_api.get_asci_symbol(idx));
        uint8_t nucleotide = kmer_api.get_num_symbol(idx);  // Works as long as our nucleotide-to-integer mapping is the same as KMC.

        kmer = (kmer << 2) | nucleotide;
    }
}


inline Kmer_u64::Kmer_u64(const Kmer_u64& rhs): kmer(rhs.kmer)
{}


inline Kmer_u64 Kmer_u64::reverse_complement() const
{
    uint64_t kmer_val = kmer;
    Kmer_u64 rev_comp;

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        rev_comp.kmer = ((rev_comp.kmer << 2) | complement_nucleotide(DNA_Base(kmer_val & 0b11)));

        kmer_val >>= 2;
    }

    return rev_comp;
}


inline void Kmer_u64::from_CKmerAPI(const CKmerAPI& kmer_api)
{
    kmer = 0;
    kmer_api.to_u64(kmer);
}


inline uint64_t Kmer_u64::to_u64() const
{
    return kmer;
}


inline bool Kmer_u64::operator<(const Kmer_u64& rhs) const
{
    return kmer < rhs.kmer;
}


inline bool Kmer_u64::operator==(const Kmer_u64& rhs) const
{
    return kmer == rhs.kmer;
}


inline bool Kmer_u64::in_forward(const Kmer_u64& kmer_hat) const
{
    return this->operator==(kmer_hat);
}


inline void Kmer_u64::roll_to_next_kmer(const char next_nucl, Kmer_u64& rev_compl)
{
    const DNA_Base mapped_nucl = map_nucleotide(next_nucl);

    kmer = ((kmer & bitmask_MSN) << 2) | mapped_nucl;
    rev_compl.kmer = (rev_compl.kmer >> 2) | (uint64_t(complement_nucleotide(mapped_nucl)) << (2 * (k - 1)));
}


inline Kmer_u64 Kmer_u64::canonical(const Kmer_u64& rev_compl) const
{
    return std::min(*this, rev_compl);
}



#endif
