
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
void createBoundaryFile(MPI_Comm comm,
                        const std::string& a_filename, int voffset, int vcount, int vtotal,
                        const std::vector<int>& a_vertexID,
                        const std::vector<P>& a_positions,
                        const std::vector<V>& a_velocity, int qoffset, int qcount, int qtotal,
                        const std::vector<int>& a_quads,
                        const std::vector<int>& a_quadID,
                        const std::vector<double>& a_pressure);

// assuming the file are opening is already laid out by createBoundaryFile
// and the VertexID and QuadIDs have not changed.
// if comm==MPI_COMM_NULL then every rank entering this function opens the file
//   ie, it not a parallel file open but a posix layer file open.
//  when "read"ing, std::vector objects get resize called on them.
void rwBoundaryFile(MPI_Comm comm, const std::string& a_filename,
                    bool a_read, double& dt,
                    int voffset, int vcount,
                    std::vector<P>& a_position,
                    std::vector<V>& a_velocity, int qoffset, int qcount,
                    std::vector<double>& a_pressure);

void rwBoundaryFile(hid_t root,
                    bool a_read, double& dt,
                    int voffset, int vcount,
                    std::vector<P>& a_position,
                    std::vector<V>& a_velocity, int qoffset, int qcount,
                    std::vector<double>& a_pressure);

// if comm == MPI_COMM_NULL then voffset and vcount are ignored and the calling
// ranks get a full copy of the data.
void readBoundaryFile(MPI_Comm comm,
                      const std::string& a_filename, double& a_dt, int voffset, int vcount,
                      std::vector<int>& a_vertexID,
                      std::vector<P>& a_positions,
                      std::vector<V>& a_velocity, int qoffset, int qcount,
                      std::vector<int>& a_quads,
                      std::vector<int>& a_quadID,
                      std::vector<double>& a_pressure);


