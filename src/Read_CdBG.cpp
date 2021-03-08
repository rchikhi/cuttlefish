
#include "Read_CdBG.hpp"


template <uint16_t k>
Read_CdBG<k>::Read_CdBG(const Build_Params& params):
    params(params)
{}


template <uint16_t k>
void Read_CdBG<k>::construct()
{
    std::cout << "\nConstructing the minimal perfect hash function (MPHF) over the vertex set.\n";
    hash_table.construct(params.vertex_db_path(), params.thread_count(), params.working_dir_path(), params.mph_file_path());


    hash_table.clear();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG)
