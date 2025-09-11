#ifdef MPI_ENABLE
    #include <mpi.h>
#endif


#include "MPI_utilities.hpp"



int getSize()
{
    int world_size;
    #ifdef MPI_ENABLE
        MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    #else
        world_size = 1;
    #endif
    return world_size;
}

int getRank()
{
    int world_rank;
    #ifdef MPI_ENABLE
        MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    #else
        world_rank = 0;
    #endif
    return world_rank;
}

// void MPI_print(std::ostream& stream, const std::string& message) ;
void MPIsafe_print(std::ostream& stream, const std::string& message) 
{
    if (getRank() == 0)
    {
        stream << message;
    }
}

void MPIsafe_exit(int exit_code)
{
	#ifdef MPI_ENABLE
		MPI_Abort(MPI_COMM_WORLD,exit_code);
    #else
        exit(exit_code);
    #endif
	
}

void MPIsafe_barrier()
{
    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
}

void MPIsafe_bcast_string(std::string& s, int root) 
{
    #ifdef MPI_ENABLE
        int len = static_cast<int>(s.size());
        MPI_Bcast(&len, 1, MPI_INT, root, MPI_COMM_WORLD);
        s.resize(len);
        if (len > 0) {
            MPI_Bcast(s.data(), len, MPI_CHAR, root, MPI_COMM_WORLD);
        }
    #endif
}