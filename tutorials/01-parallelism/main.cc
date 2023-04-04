#include "spade.h"

int main(int argc, char** argv)
{
    using real_t = double;

    //Initialize MPI using SPADE's MPI wrapper
    spade::parallel::mpi_t group(&argc, &argv);
    
    //Root rank does some printing
    if (group.isroot()) print("Initialized parallel group.");
    if (group.isroot()) print("Size:", group.size());

    //Each rank will generate some data:
    spade::utils::random_seed(group.rank());
    std::vector<real_t> data(group.rank()+2, 0.0);
    for (auto& e: data) e = spade::utils::unitary_random(); //Generates a random number between 0 and 1

    if (group.isroot()) print("Report from each rank (with some data)");

    //SPADE currently provides a range iterator for convenience, but note that it is not considered best practice!
    group.sync();
    for (auto qrank: range(0, group.size()))
    {
        if (qrank == group.rank())
        {
            print("Rank:", group.rank());
            print("Data:");
            for (const auto& e: data)
            {
                print(e);
            }
            print();
        }

        //Wrapper for MPI_Barrier
        group.sync();
    }
    group.sync();

    //We can sum up values using the group:
    auto total_data_size = group.sum(data.size());
    if (group.isroot()) print("Total data size:", total_data_size);

    //We can take the total data sum as well:
    auto total_data_sum  = group.sum(data);
    if (group.isroot()) print("Total data sum:", total_data_sum);

    //See ${SPADE}/src/parallel.h for more details.

    return 0;
}
