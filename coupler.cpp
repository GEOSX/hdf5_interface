#include "coupler.H"

namespace internal
{

/*!
 * \brief Create a dataset under the given group with the given data.
 *
 * \param [in] group the group where the dataset will live.
 * \param [in] name the name of the dataset.
 * \param [in] type the datatype of the dataset.
 * \param [in] global_offset the offset at which this rank is to write to
 *  the dataset.
 * \param [in] global_size the total size of the dataset.
 * \param [in] local_size the size of the chunk to be written by this rank.
 * \param [in] data the data to be written.
 *
 * \note all indexing is done in units of TYPE. So if TYPE is a vector of three
 *  doubles then a LOCAL_SIZE of 4 means that you write a total of 12 doubles.
 */
void createDataset(hid_t group, const char* name, hid_t type, 
                   hsize_t global_offset, hsize_t global_size,
                   hsize_t local_size, const void* data)
{
  hid_t dataspace = H5Screate_simple(1, &global_size, NULL);
  hid_t dataset = H5Dcreate(group, name, type, dataspace, H5P_DEFAULT, 
                            H5P_DEFAULT, H5P_DEFAULT);

  hid_t mem_dataspace = H5Screate_simple(1, &local_size, NULL);
  hid_t hyperslab = dataspace;
  H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &global_offset, NULL, 
                      &local_size, NULL);
  H5Dwrite(dataset, type, mem_dataspace, hyperslab, H5P_DEFAULT, data);
  
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(mem_dataspace);
}

/*!
 * \brief Create a dataset under the given group with a subset of the given data
 *  specified by the indices arra.
 *
 * \param [in] group the group where the dataset will live.
 * \param [in] name the name of the dataset.
 * \param [in] type the datatype of the dataset.
 * \param [in] global_offset the offset at which this rank is to write to the dataset.
 * \param [in] global_size the total size of the dataset.
 * \param [in] local_size the size of the chunk to be written by this rank.
 * \param [in] buffer_size the total length of the DATA array.
 * \param [in] indices the indices of the values in the DATA array to write out.
 * \param [in] data the data to be written.
 *
 * \note all indexing is done in units of TYPE. So if TYPE is a vector of three
 *  doubles then a LOCAL_SIZE of 4 means that you write a total of 12 doubles.
 */
void createDataset(hid_t group, const char* name, hid_t type, 
                   hsize_t global_offset, hsize_t global_size,
                   hsize_t local_size, hsize_t buffer_size, 
                   const hsize_t* indices,
                   const void* data)
{
  hid_t dataspace = H5Screate_simple(1, &global_size, NULL);
  hid_t dataset = H5Dcreate(group, name, type, dataspace, H5P_DEFAULT, 
                            H5P_DEFAULT, H5P_DEFAULT);

  hid_t mem_dataspace = H5Screate_simple(1, &buffer_size, NULL);
  H5Sselect_elements(mem_dataspace, H5S_SELECT_SET, local_size, indices);
  
  hid_t hyperslab = dataspace;
  H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &global_offset, NULL, 
                      &local_size, NULL);
  H5Dwrite(dataset, type, mem_dataspace, hyperslab, H5P_DEFAULT, data);
  
  H5Sclose(mem_dataspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

/*!
 * \brief Create and write field data.
 *
 * \param [in] group the group under which to create the datasets.
 * \param [in] fm the fields to write out.
 * \param [in] global_offset the offset at which this rank is to write.
 * \param [in] global_size the total size of the dataset.
 * \param [in] local_size the size of the chunk to be written by this rank.
 * \param [in] buffer_size the total length of the DATA array.
 * \param [in] indices the indices of the values in the DATA array to write out.
 */
void createFields(hid_t group, const FieldMap& fm, hsize_t global_offset, 
                  hsize_t global_size, hsize_t local_size, 
                  hsize_t buffer_size, const hsize_t* indices)
{
  for (FieldMap::const_iterator it = fm.begin(); 
       it != fm.end(); ++it)
  {
    /* Extract the field data and metadata. */
    const char* name = it->first.c_str();
    const hid_t base_type = std::get<0>(it->second);
    const hsize_t n_components = std::get<1>(it->second);
    const void* data = std::get<2>(it->second);

    /* Create the appropriate type and dataset. */
    hid_t field_type = H5Tarray_create(base_type, 1, &n_components);
    internal::createDataset(group, name, field_type, global_offset, 
                            global_size, local_size, 
                            buffer_size, indices, data);
    H5Tclose(field_type);
  }
}

/*!
 * \brief Read a dataset into the given array.
 *
 * \param [in] group the group where the dataset lives.
 * \param [in] name the name of the dataset.
 * \param [in] global_offset the offset at which this rank is to read from the 
 *  dataset.
 * \param [in] n_elems the number of elements to read.
 * \param [out] data the array to be read into.
 *
 * \note all indexing is done in units of the dataset type. So if the type is a
 *  vector of three doubles then an N_ELEMS of 4 means that you read a total of 
 *  12 doubles.
 */
void readDataset(hid_t group, const char* name, hsize_t global_offset, 
                 hsize_t n_elems, void* data)
{
  /* Open the dataset and get the type. */
  hid_t dataset = H5Dopen(group, name, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);
  hid_t type = H5Dget_type(dataset);

  /* Check that the dataset is 1D. */
  int ndims = H5Sget_simple_extent_ndims(dataspace);
  if (ndims != 1)
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* Create the memory and file dataspaces and do the read. */
  hid_t mem_dataspace = H5Screate_simple(1, &n_elems, NULL);
  hid_t hyperslab = dataspace;
  H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &global_offset, 
                      NULL, &n_elems, NULL);

  H5Dread(dataset, type, mem_dataspace, hyperslab, H5P_DEFAULT, data);

  H5Sclose(dataspace);
  H5Sclose(mem_dataspace);
  H5Tclose(type);
  H5Dclose(dataset);
}

/*!
 * \brief Read a dataset into the given positions of the given array.
 *
 * \param [in] group the group where the dataset lives.
 * \param [in] name the name of the dataset.
 * \param [in] global_offset the offset at which this rank is to read from the 
 *  dataset.
 * \param [in] n_elems the number of elements to read.
 * \param [in] buffer_size the size of the DATA array.
 * \param [in] indices the indices of the values in the DATA array to read to.
 * \param [out] data the array to be read into.
 *
 * \note all indexing is done in units of the dataset type. So if the type is a
 *  vector of three doubles then an N_ELEMS of 4 means that you read a total of 
 *  12 doubles.
 */
void readDataset(hid_t group, const char* name, hsize_t global_offset, 
                 hsize_t n_elems, hsize_t buffer_size, const hsize_t* indices,
                 void* data)
{
  /* Open the dataset and get the type. */
  hid_t dataset = H5Dopen(group, name, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);
  hid_t type = H5Dget_type(dataset);

  /* Check that the dataset is 1D. */
  int ndims = H5Sget_simple_extent_ndims(dataspace);
  if (ndims != 1)
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* Create the memory and file dataspaces and do the read. */
  hid_t mem_dataspace = H5Screate_simple(1, &buffer_size, NULL);
  H5Sselect_elements(mem_dataspace, H5S_SELECT_SET, n_elems, indices);

  hid_t hyperslab = dataspace;
  H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, &global_offset, 
                      NULL, &n_elems, NULL);

  H5Dread(dataset, type, mem_dataspace, hyperslab, H5P_DEFAULT, data);

  H5Sclose(dataspace);
  H5Sclose(mem_dataspace);
  H5Tclose(type);
  H5Dclose(dataset);
}

/*!
 * \brief Read field data.
 *
 * \param [in] group the group under which the fields live.
 * \param [in/out] fm the fields to write to.
 * \param [in] global_offset the offset at which this rank is to write.
 * \param [in] n_elems the number of elements to read in.
 * \param [in] buffer_size the size of the field data arrays.
 * \param [in] indices the indices of the values in the field arrays to write to.
 */
void readFields(hid_t group, FieldMap& fm, hsize_t global_offset, 
                hsize_t n_elems, hsize_t buffer_size, const hsize_t* indices)
{
  for (FieldMap::iterator it = fm.begin(); 
       it != fm.end(); ++it)
  {
    /* Extract the field data and metadata. */
    const char* name = it->first.c_str();
    void* data = std::get<2>(it->second);

    /* Read the dataset. */
    internal::readDataset(group, name, global_offset, n_elems, buffer_size, 
                          indices, data);
  }
}

/*!
 * \brief Get the offset at which to write to.
 *
 * \param [in] comm the communicator used to write the boundary file.
 * \param [in] local_count the number of elements this rank will write out.
 * \param [in] offset the offset at which this rank is to write.
 * \param [in] global_count the total number of elements in the dataset.
 */
void boundaryFileOffsets(MPI_Comm comm, int local_count, int& offset, 
                         int& global_count)
{
  int N;
  MPI_Comm_size(comm, &N);
  MPI_Scan(&local_count, &offset, 1, MPI_INT, MPI_SUM, comm );
  global_count = offset;
  MPI_Bcast( &global_count, 1, MPI_INT, N-1, comm); //last rank has total sum
  offset -= local_count;
}

} /* namespace internal */


//------------------------------------------------------------------------------
void writeBoundaryFile(MPI_Comm comm, const char* filename, int& face_offset,
                       int& n_faces_to_write, int& node_offset, 
                       int& n_nodes_to_write, double dt, int n_faces,
                       int n_nodes, const int* faces, const bool* on_boundary,
                       const FieldMap& face_fields, const FieldMap& node_fields)
{
  /* Create the file and open the root group. */
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  if (comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);
  }

  hid_t m_fileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, 
                             file_access);
  H5Pclose(file_access);
  hid_t root = H5Gopen(m_fileID, "/", H5P_DEFAULT);

  /* Write dt. */
  hid_t aid  = H5Screate(H5S_SCALAR);
  hid_t attr1 = H5Acreate(root,"dt",H5T_NATIVE_DOUBLE, aid, H5P_DEFAULT, 
                          H5P_DEFAULT);
  H5Awrite(attr1, H5T_NATIVE_DOUBLE, &dt);
  H5Aclose(attr1);
  H5Sclose(aid);

  /* Calculate the number of faces to write. */
  n_faces_to_write = 0;
  for (int i = 0; i < n_faces; ++i)
  {
    n_faces_to_write += on_boundary[i];
  }

  /* Get the face offset at which to write and the total number of faces to write. */
  face_offset = 0;
  int total_n_faces_to_write = n_faces_to_write;
  internal::boundaryFileOffsets(comm, n_faces_to_write, face_offset, 
                                total_n_faces_to_write);

  /* Collect the face local IDs to write. */
  hsize_t* face_IDs = new hsize_t[n_faces_to_write];
  int cur_n_faces = 0;
  for (int i = 0; i < n_faces; ++i)
  {
    if (on_boundary[i])
    {
      face_IDs[cur_n_faces] = i;
      cur_n_faces++;
    }
  }

  /* Create the map of local node IDs to file node IDs and the reverse map. */
  std::map<int, int> nodeLocalToFileID;
  for (int i = 0; i < n_faces_to_write; ++i)
  {
    const int cur_face = face_IDs[i];
    for (int j = 0; j < 4; ++j)
    {
      int cur_nodeID = faces[4 * cur_face + j];
      nodeLocalToFileID[cur_nodeID] = -1;
    }
  }

  int fileNodeID = 0;
  n_nodes_to_write = nodeLocalToFileID.size();
  hsize_t* nodeFileToLocalID = new hsize_t[n_nodes_to_write];
  for (std::map<int, int>::iterator it = nodeLocalToFileID.begin();
       it != nodeLocalToFileID.end(); ++it)
  {
    nodeFileToLocalID[fileNodeID] = it->first;
    it->second = fileNodeID;
    fileNodeID++;
  }

  /* Get the node offset at which to write and the total number of nodes to write. */
  node_offset = 0;
  int total_n_nodes_to_write = n_nodes_to_write;
  internal::boundaryFileOffsets(comm, n_nodes_to_write, node_offset, 
                                total_n_nodes_to_write);

  /* Write out the node local IDs. */
  internal::createDataset(root, "NodeIDs", H5T_NATIVE_HSIZE, node_offset,
                          total_n_nodes_to_write, n_nodes_to_write, 
                          nodeFileToLocalID);

  /* Write out the node field data. */
  internal::createFields(root, node_fields, node_offset, 
                         total_n_nodes_to_write, n_nodes_to_write,
                         n_nodes, nodeFileToLocalID);

  delete[] nodeFileToLocalID;
  nodeFileToLocalID = nullptr;

  /* Redo the face connectivity with the file nodal numbering. */
  int* faces_to_write = new int[4 * n_faces_to_write];
  for (int i = 0; i < n_faces_to_write; ++i)
  {
    const hsize_t cur_face = face_IDs[i];
    for (int j = 0; j < 4; ++j)
    {
      const int cur_node = faces[4 * cur_face + j];
      faces_to_write[4 * i + j] = nodeLocalToFileID[cur_node];
    }
  }

  /* Write out the faces. */
  const hsize_t verts_per_quad = 4;
  hid_t quad_type = H5Tarray_create(H5T_NATIVE_INT, 1, &verts_per_quad);
  internal::createDataset(root, "Quads", quad_type, face_offset,
                          total_n_faces_to_write, n_faces_to_write, 
                          faces_to_write);
  H5Tclose(quad_type);

  delete faces_to_write;
  faces_to_write = nullptr;

  /* Write out the face local IDs. */
  internal::createDataset(root, "QuadIDs", H5T_NATIVE_HSIZE, face_offset, 
                          total_n_faces_to_write, n_faces_to_write, face_IDs);

  /* Write out the face field data. */
  internal::createFields(root, face_fields, face_offset, 
                         total_n_faces_to_write, n_faces_to_write, 
                         n_faces, face_IDs);

  delete[] face_IDs;
  face_IDs = nullptr;  

  H5Gclose(root);
  H5Fclose(m_fileID);
}

//------------------------------------------------------------------------------
void readBoundaryFile(MPI_Comm comm, const char* filename, int face_offset,
                      int n_faces_to_read, int n_faces, int node_offset, 
                      int n_nodes_to_read, int n_nodes, FieldMap& face_fields,
                      FieldMap& node_fields)
{
  /* Open the file and the root group. */
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  if(comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);
  }

  hid_t m_fileID = H5Fopen(filename, H5F_ACC_RDONLY, file_access);
  H5Pclose(file_access);
  hid_t root = H5Gopen(m_fileID, "/", H5P_DEFAULT);

  /* Get the localIDs of the faces that were written. */
  hsize_t* face_IDs = new hsize_t[n_faces_to_read];
  internal::readDataset(root, "QuadIDs", face_offset, n_faces_to_read, face_IDs);

  /* Read in the face fields. */
  internal::readFields(root, face_fields, face_offset, n_faces_to_read, n_faces, 
                       face_IDs);
  delete[] face_IDs;

  /* Get the localIDs of the nodes that were written. */
  hsize_t* node_IDs = new hsize_t[n_nodes_to_read];
  internal::readDataset(root, "NodeIDs", node_offset, n_nodes_to_read, node_IDs);

  /* Read in the node fields. */
  internal::readFields(root, node_fields, node_offset, n_nodes_to_read, n_nodes, 
                       node_IDs);
  delete[] node_IDs;

  H5Gclose(root);
  H5Fclose(m_fileID);
}