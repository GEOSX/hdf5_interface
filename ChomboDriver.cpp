
#include "coupler.H"


int stop(int, void*){ abort(); return 0;}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv );
  int rank, N;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &N);

  H5Eset_auto (H5E_DEFAULT,stop, NULL); //for debugging, tell HDF5 to abort on error

  int VCOUNT = 0;
  int QCOUNT = 0;
  int* VID = nullptr;   
  P* x = nullptr;
  V* v = nullptr;
  double* pressure = nullptr;
  int* QID = nullptr;
  int* quads = nullptr;
 
  double dt;

  if(rank == 0)
  {
    readBoundaryFile(MPI_COMM_NULL, "GEOSboundary.hdf5", dt, 0, VCOUNT, VID, x,
                     v, 0, QCOUNT, quads, QID, pressure);
  }
  
  for(int i = 0; i < QCOUNT; ++i) 
  {
    pressure[i] += 0.055;
  }
  
  // int d[2] = {(int)VID.size(), (int)QID.size()};
  // MPI_Bcast(d, 2, MPI_INT, 0, MPI_COMM_WORLD);
  // VID.resize(d[0]);
  // QID.resize(d[1]);
  // x.resize(d[0]);
  // v.resize(d[0]);
  // pressure.resize(d[1]);
  // MPI_Bcast(pressure.data(), d[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  V* empty = nullptr;
  P* emptyP = nullptr;

  dt = 1.155;
  // rwBoundaryFile(MPI_COMM_WORLD, "GEOSboundary.hdf5.tmp2", false, //false="write"
  //                dt, 0, 0, emptyP, empty, quadOffset, QCOUNT, pressure);
  if (rank==0)
  {
    writeBoundaryFile(MPI_COMM_NULL, "GEOSboundary.hdf5", dt, 0, VCOUNT, emptyP,
                      empty, 0, QCOUNT, pressure);
  }

  
  
  MPI_Finalize();
  
  return 0;
}
