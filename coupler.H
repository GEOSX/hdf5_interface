
#include <mpi.h>
#include <hdf5.h>
#include <vector>
#include <string>
#include <cstdio>

struct P {
  double X, Y, Z;
};

struct V {
  double X, Y, Z;
};

void boundaryFileOffsets(MPI_Comm comm, int localCount, int& offset, int& count);

// vectors here are local portion of total boundary arrays
void createBoundaryFile(MPI_Comm comm, const std::string& filename, int voffset,
                        int vcount, int vtotal, const int* vertexID, 
                        const P* positions, const V* velocity, int qoffset,
                        int qcount, int qtotal, const int* quads, 
                        const int* quadID, const double* pressure);

// assuming the file are opening is already laid out by createBoundaryFile
// and the VertexID and QuadIDs have not changed.
// if comm==MPI_COMM_NULL then every rank entering this function opens the file
//   ie, it not a parallel file open but a posix layer file open.
//  when "read"ing, std::vector objects get resize called on them.
void readBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt,
                      int voffset, int vcount, P*& positions, V*& velocities, 
                      int qoffset, int qcount, double*& pressures);

void readBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt, 
                      int voffset, int& vcount, int*& vertexIDs, P*& positions,
                      V*& velocities, int qoffset, int& qcount, int*& quads,
                      int*& quadIDs, double*& pressures);

void writeBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt,
                       int voffset, int vcount, const P* positions, 
                       const V* velocities, int qoffset, int qcount,
                       const double* pressures);

