#include "coupler.H"


namespace internal
{

void readBoundaryFile(hid_t root, double& dt, int voffset, int vcount, 
                      P*& positions, V*& velocities, int qoffset, int qcount, 
                      double*& pressures)
{
  hid_t attr1= H5Aopen(root, "dt", H5P_DEFAULT);
  hid_t aid = H5Aget_space(attr1);
  H5Aread(attr1, H5T_NATIVE_DOUBLE, &dt);

  hsize_t dims[1];
  hsize_t offset[1];
  hsize_t count[1];
  
  hid_t pressures_d = H5Dopen(root, "Pressure", H5P_DEFAULT);
  hid_t pressures_s = H5Dget_space(pressures_d);
  dims[0] = qcount;
  hid_t mdataspace = H5Screate_simple(1, dims, NULL);
  hid_t hyperslab = pressures_s;
  if (qcount > 0)
  {
    offset[0] = qoffset;
    count[0] = qcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset , NULL, count, NULL);
  }
  else
  {
    H5Sselect_none(hyperslab);
  }

  delete[] pressures;
  pressures = new double[qcount];
  H5Dread(pressures_d, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
          pressures);

  hid_t p_set = H5Dopen(root, "Position", H5P_DEFAULT);
  hid_t p_space = H5Dget_space(p_set);
  hyperslab = p_space;

  dims[0] = 3 * vcount;
  mdataspace = H5Screate_simple(1, dims, NULL);

  if (vcount > 0)
  {
    offset[0] = 3 * voffset;
    count[0] = 3 * vcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset , NULL, count, NULL);
  }
  else
  {
    H5Sselect_none(hyperslab);
  }

  delete positions;
  delete[] velocities;
  positions = new P[vcount];
  velocities = new V[vcount];
  H5Dread(p_set, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
          positions);
  hid_t v_set = H5Dopen(root, "Velocity", H5P_DEFAULT);
  H5Dread(v_set,  H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
          velocities);

  H5Sclose(mdataspace);
  H5Sclose(pressures_s);
  H5Sclose(aid);
  H5Dclose(p_set);
  H5Sclose(p_space);
  H5Dclose(v_set);
  H5Dclose(pressures_d);
  H5Aclose(attr1);
}


void writeBoundaryFile(hid_t root, double& dt, int voffset, int vcount, 
                       const P* positions, const V* velocities, int qoffset, 
                       int qcount, const double* pressures)
{
  hid_t attr1= H5Aopen(root, "dt", H5P_DEFAULT);
  hid_t aid = H5Aget_space(attr1);
  H5Awrite(attr1, H5T_NATIVE_DOUBLE, &dt);

  hsize_t dims[1];
  hsize_t offset[1];
  hsize_t count[1];
  
  hid_t pressures_d = H5Dopen(root, "Pressure", H5P_DEFAULT);
  hid_t pressures_s = H5Dget_space(pressures_d);
  dims[0] = qcount;
  hid_t mdataspace = H5Screate_simple(1, dims, NULL);
  hid_t hyperslab = pressures_s;
  if (qcount > 0)
  {
    offset[0] = qoffset;
    count[0] = qcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset , NULL, count, NULL);
  }
  else
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(pressures_d, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT,
           pressures);

  hid_t p_set = H5Dopen(root, "Position", H5P_DEFAULT);
  hid_t p_space = H5Dget_space(p_set);
  hyperslab = p_space;

  dims[0] = 3 * vcount;
  mdataspace = H5Screate_simple(1, dims, NULL);

  if (vcount > 0)
  {
    offset[0] = 3 * voffset;
    count[0] = 3 * vcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset , NULL, count, NULL);
  }
  else
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(p_set, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
           positions);
  hid_t v_set = H5Dopen(root, "Velocity", H5P_DEFAULT);
  H5Dwrite(v_set,  H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
           velocities);

  H5Sclose(mdataspace);
  H5Sclose(pressures_s);
  H5Sclose(aid);
  H5Dclose(p_set);
  H5Sclose(p_space);
  H5Dclose(v_set);
  H5Dclose(pressures_d);
  H5Aclose(attr1);
}

} /* namespace internal */


void boundaryFileOffsets(MPI_Comm comm, int localCount, int& offset, int& count)
{
  int N;
  MPI_Comm_size(comm, &N);
  MPI_Scan (&localCount, &offset, 1, MPI_INT, MPI_SUM, comm );
  count = offset;
  MPI_Bcast( &count, 1, MPI_INT, N-1, comm); //last rank has total sum
  offset -= localCount;
}

// vectors here are local portion of total boundary arrays
void createBoundaryFile(MPI_Comm comm, const std::string& filename, int voffset,
                        int vcount, int vtotal, const int* vertexIDs, 
                        const P* positions, const V* velocity, int qoffset,
                        int qcount, int qtotal, const int* quads, 
                        const int* quadIDs, const double* pressures)
{
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  if (comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);
  }
  
  hid_t m_fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, 
                             file_access);
  H5Pclose(file_access);

  hid_t root = H5Gopen(m_fileID, "/", H5P_DEFAULT);

  double dt=0.0;
  hid_t aid  = H5Screate(H5S_SCALAR);
  hid_t attr1 = H5Acreate(root,"dt",H5T_NATIVE_DOUBLE, aid, H5P_DEFAULT, 
                          H5P_DEFAULT);
  H5Awrite(attr1, H5T_NATIVE_DOUBLE, &dt);
  H5Sclose(aid);
  
  hsize_t dims[1];
  dims[0] = vtotal;
  hid_t fdataspace = H5Screate_simple(1, dims, NULL);
  hid_t id_set = H5Dcreate(root, "VertexID", H5T_NATIVE_INT, fdataspace, 
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t hyperslab = fdataspace;

  dims[0] = vcount;
  hid_t mdataspace = H5Screate_simple(1, dims, NULL);

  hsize_t offset[1];
  hsize_t count[1];
  if (vcount > 0)
  {
    offset[0] = voffset;
    count[0] = vcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offset , NULL, count, NULL);
  } 
  else 
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(id_set, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, 
           vertexIDs);
  
  H5Sclose(mdataspace);
  H5Sclose(fdataspace);
  dims[0] = 3 * vtotal;
  fdataspace = H5Screate_simple(1, dims, NULL);
  hid_t p_set = H5Dcreate(root, "Position", H5T_NATIVE_DOUBLE, fdataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hyperslab = fdataspace;

  dims[0]= 3 * vcount;
  mdataspace = H5Screate_simple(1, dims, NULL);

  if (vcount > 0)
  {
    offset[0]=3*voffset;
    count[0]=3*vcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET,
                        offset , NULL,
                        count, NULL);
  } 
  else 
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(p_set, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
           positions);
  hid_t v_set = H5Dcreate(root, "Velocity", H5T_NATIVE_DOUBLE, fdataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(v_set,  H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
           velocity);

  
  dims[0] = qtotal * 4;
  fdataspace = H5Screate_simple(1, dims, NULL);
  hid_t q_set = H5Dcreate(root, "Quads", H5T_NATIVE_INT, fdataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hyperslab = fdataspace;

  dims[0] = qcount * 4;
  mdataspace = H5Screate_simple(1, dims, NULL);

  if (qcount > 0)
  {
    offset[0]=4*qoffset;
    count[0]=4*qcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET,
                        offset , NULL,
                        count, NULL);
  } 
  else
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(q_set, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, quads);

  dims[0]= qtotal;
  fdataspace = H5Screate_simple(1, dims, NULL);
  hid_t pres = H5Dcreate(root, "Pressure", H5T_NATIVE_DOUBLE, fdataspace, 
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t qid  = H5Dcreate(root, "QuadID", H5T_NATIVE_INT, fdataspace, 
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hyperslab = fdataspace;

  dims[0] = qcount;
  mdataspace = H5Screate_simple(1, dims, NULL);

  if (qcount > 0)
  {
    offset[0]=qoffset;
    count[0]=qcount;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET,
                        offset , NULL,
                        count, NULL);
  }
  else
  {
    H5Sselect_none(hyperslab);
  }

  H5Dwrite(pres, H5T_NATIVE_DOUBLE, mdataspace, hyperslab, H5P_DEFAULT, 
           pressures);

  H5Dwrite(qid, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, quadIDs);
 
  H5Sclose(mdataspace);
  H5Sclose(fdataspace);
  H5Dclose(p_set);
  H5Dclose(id_set);
  H5Dclose(v_set);
  H5Dclose(q_set);
  H5Dclose(pres);
  H5Dclose(qid);
  H5Aclose(attr1);
  H5Gclose(root);
  H5Fclose(m_fileID);
}

// assuming the file are opening is already laid out by createBoundaryFile
// and the VertexID and QuadIDs have not changed.
// if comm==MPI_COMM_NULL then every rank entering this function opens the file
//   ie, it not a parallel file open but a posix layer file open.
//  when "read"ing, std::vector objects get resize called on them.
void readBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt,
                      int voffset, int vcount, P*& positions, V*& velocities, 
                      int qoffset, int qcount, double*& pressures)
{
  hid_t file_access = H5Pcreate (H5P_FILE_ACCESS);
  if(comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);
  }

  hid_t m_fileID = H5Fopen(filename.c_str(), H5F_ACC_RDWR, file_access);
  H5Pclose(file_access);
  hid_t root = H5Gopen(m_fileID, "/", H5P_DEFAULT);
  internal::readBoundaryFile(root, dt, voffset, vcount, positions, velocities, 
                             qoffset, qcount, pressures);
  H5Gclose(root);
  H5Fclose(m_fileID);
}

void readBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt, 
                      int voffset, int& vcount, int*& vertexIDs, P*& positions,
                      V*& velocities, int qoffset, int& qcount, int*& quads,
                      int*& quadIDs, double*& pressures)
{
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  if(comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access,  comm, MPI_INFO_NULL);
  }

  hid_t m_fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, file_access);
  H5Pclose(file_access);
  hid_t root = H5Gopen(m_fileID, "/",H5P_DEFAULT);

  hid_t vid_d = H5Dopen(root, "VertexID", H5P_DEFAULT);
  hid_t vid_s = H5Dget_space(vid_d);
  hsize_t dims[1];
  hsize_t offsets[1];
  H5Sget_simple_extent_dims(vid_s, dims, NULL);
  hid_t mdataspace, hyperslab;
  if(comm == MPI_COMM_NULL)
  {
    mdataspace = H5Screate_simple(1, dims, NULL);
    hyperslab = mdataspace;
    delete[] vertexIDs;
    vertexIDs = new int[dims[0]];
    voffset = 0;
    vcount = dims[0];
  }
  else
  {
    dims[0] = vcount;
    delete[] vertexIDs;
    vertexIDs = new int[vcount];
    mdataspace = H5Screate_simple(1, dims, NULL);
    hyperslab = mdataspace;
    if (vcount > 0)
    {
      offsets[0] = voffset;
      H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET,
                          offsets , NULL,
                          dims, NULL);
    }
    else
    {
      H5Sselect_none(hyperslab);
    }
  }

  H5Dread(vid_d, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, vertexIDs);
  H5Sclose(mdataspace);
  
  hid_t qid_d = H5Dopen(root, "QuadID", H5P_DEFAULT);
  hid_t qid_s = H5Dget_space(qid_d);
  H5Sget_simple_extent_dims(qid_s, dims, NULL);
  if (comm == MPI_COMM_NULL)
  {
    mdataspace = H5Screate_simple(1, dims, NULL);
    hyperslab = mdataspace;
    delete[] quadIDs;
    quadIDs = new int[dims[0]];
    qoffset = 0;
    qcount = dims[0];
  } 
  else
  {
    dims[0] = qcount;
    delete[] quadIDs;
    quadIDs  = new int[qcount];
    mdataspace = H5Screate_simple(1, dims, NULL);
    hyperslab=mdataspace;
    if(qcount > 0)
    {
      offsets[0] = qoffset;
      H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offsets , NULL, dims, NULL);
    } 
    else 
    {
      H5Sselect_none(hyperslab);
    }
  }

  H5Dread(qid_d, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, quadIDs);
  H5Sclose(mdataspace);

  delete[] quads;
  quads = new int[4 * qcount];
  hid_t quad_d = H5Dopen(root, "Quads", H5P_DEFAULT);
  hid_t quad_s = H5Dget_space(quad_d);
  dims[0] = 4 * qcount;
  mdataspace = H5Screate_simple(1, dims, NULL);
  hyperslab = mdataspace;
  if (qcount > 0)
  {
    offsets[0] = qoffset;
    H5Sselect_hyperslab(hyperslab, H5S_SELECT_SET, offsets , NULL, dims, NULL);
  }
  else 
  {
    H5Sselect_none(hyperslab);
  }
  H5Dread(quad_d, H5T_NATIVE_INT, mdataspace, hyperslab, H5P_DEFAULT, quads);

  internal::readBoundaryFile(root, dt, voffset, vcount, positions, velocities, 
                            qoffset, qcount, pressures);

  H5Sclose(mdataspace);
  H5Sclose(vid_s);
  H5Sclose(qid_s);
  H5Sclose(quad_s);
  H5Dclose(vid_d);
  H5Dclose(qid_d);
  H5Dclose(quad_d);
  H5Gclose(root);
  H5Fclose(m_fileID);
}


// assuming the file are opening is already laid out by createBoundaryFile
// and the VertexID and QuadIDs have not changed.
// if comm==MPI_COMM_NULL then every rank entering this function opens the file
//   ie, it not a parallel file open but a posix layer file open.
void writeBoundaryFile(MPI_Comm comm, const std::string& filename, double& dt,
                       int voffset, int vcount, const P* positions, 
                       const V* velocities, int qoffset, int qcount,
                       const double* pressures)
{
  hid_t file_access = H5Pcreate (H5P_FILE_ACCESS);
  if(comm != MPI_COMM_NULL)
  {
    H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);
  }

  hid_t m_fileID = H5Fopen(filename.c_str(), H5F_ACC_RDWR, file_access);
  H5Pclose(file_access);
  hid_t root = H5Gopen(m_fileID, "/", H5P_DEFAULT);
  internal::writeBoundaryFile(root, dt, voffset, vcount, positions, velocities, 
                               qoffset, qcount, pressures);
  H5Gclose(root);
  H5Fclose(m_fileID);
}


