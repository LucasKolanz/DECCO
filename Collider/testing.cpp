//TODO: Restarting a relax job double copies the constant file into the RELAXconstants file, doubling its size and throwing off
//      Porosity_FD.py. The double copy happens for RELAXenergy and RELAXsimData as well, but Porosity_FD only looks at the last line, 
//      which is still good data. Should probably just delete all RELAX files that already exist for a RELAX restart. 



// This file was originally written for SpaceLab/DECCO take care of initialization and calling the correct functions from ball_group.hpp to 
// carry out Soft Sphere, Discrete Element simulations 

// Authors: Job Guidos and Lucas Kolanz


#include <omp.h>


#include "ball_group.hpp"
#include "../timing/timing.hpp"
#include "../utilities/Utils.hpp"

#ifdef MPI_ENABLE
    #include <mpi.h>
#endif


//How many lines will be saved to a buffer before being saved
extern const int bufferlines;


// // Prototypes
// void
// safetyChecks(Ball_group &O);
// void 
// runAggregation(std::string path, int num_balls);
// void 
// runRelax(std::string path);
// void 
// runCollider(std::string path);
// timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(int argc, char* argv[])
{
    std::string path = "/home/kolanzl/Desktop/Visualize/V19/";
    Ball_group O = Ball_group(path,true);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

