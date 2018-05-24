#include <mpi.h>
#include <hdf5.h>
#include <vector>
#include <string>
#include <cstdio>
#include <tuple>
#include <map>

struct vec3
{
  double e0;
  double e1;
  double e2;
};

/* Map from field name to the HDF5 data type, the number of components per
 * object and a pointer to the data. */
using FieldMap = std::map<std::string, std::tuple<hid_t, int, void*>>;

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
void writeBoundaryFile(MPI_Comm comm, const char* filename, int& face_offset,
                       int& n_faces_to_write, int& node_offset, 
                       int& n_nodes_to_write, double dt, int n_faces,
                       int n_nodes, const int* faces, const bool* on_boundary,
                       const FieldMap& face_fields, const FieldMap& node_fields);


void readBoundaryFile(MPI_Comm comm, const char* filename, int face_offset,
                      int n_faces_to_read, int n_faces, int node_offset, 
                      int n_nodes_to_read, int n_nodes, FieldMap& face_fields,
                      FieldMap& node_fields);

