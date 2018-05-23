#include <mpi.h>
#include <hdf5.h>
#include <vector>
#include <string>
#include <cstdio>
#include <tuple>
#include <unordered_map>

/* Map from field name to the HDF5 data type, the number of components per
 * object and a pointer to the data. */
using FieldMap = std::unordered_map<std::string, std::tuple<hid_t, int, void*>>;

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

/*!
 * \brief write out a boundary file with the given data.
 *
 * \param [in] comm the communicator used in writing the file.
 * \param [in] filename the name of the file to write out to.
 * \param [in] dt the current time step.
 * \param [in] n_nodes the number of nodes in the entire mesh.
 * \param [in] n_faces the number of faces in the entire mesh.
 * \param [in] faces connectivity of the faces which are assumed to be quads.
 * \param [in] on_boundary true iff the respective face is to be written out.
 * \param [in] face_fields map from face field names to fields.
 * \param [in] node_fields map from node field names to fields.
 */
void writeBoundaryFile(MPI_Comm comm, const std::string& filename, double dt,
                       int n_nodes, int n_faces, const int* faces, 
                       const bool* on_boundary, const FieldMap& face_fields,
                       const FieldMap& node_fields)

