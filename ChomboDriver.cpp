
#include "coupler.H"


int stop(int, void*){ abort(); return 0;}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv );
  int rank, N;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &N);

  H5Eset_auto (H5E_DEFAULT,stop, NULL); //for debugging, tell HDF5 to abort on error

  std::vector<int> VID;   
  std::vector<P> x;
  std::vector<V> v;
  std::vector<double> pressure;
  std::vector<int> QID;
  std::vector<int> quads;
 
  double dt;

  if(rank == 0)
    readBoundaryFile(MPI_COMM_NULL, "GEOSboundary.hdf5", dt, 0, 0, VID,
                     x, v, 0, 0, quads, QID, pressure);
  
  for(auto p=pressure.begin(); p!= pressure.end(); p++) *p+=0.055;
  
  // int d[2] = {(int)VID.size(), (int)QID.size()};
  // MPI_Bcast(d, 2, MPI_INT, 0, MPI_COMM_WORLD);
  // VID.resize(d[0]);
  // QID.resize(d[1]);
  // x.resize(d[0]);
  // v.resize(d[0]);
  // pressure.resize(d[1]);
  // MPI_Bcast(pressure.data(), d[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  std::vector<V> empty;
  std::vector<P> emptyP;

  dt = 1.155;
  // rwBoundaryFile(MPI_COMM_WORLD, "GEOSboundary.hdf5.tmp2", false, //false="write"
  //                dt, 0, 0, emptyP, empty, quadOffset, QCOUNT, pressure);
  if(rank==0)
    {
      rwBoundaryFile(MPI_COMM_NULL, "GEOSboundary.hdf5", false, //false="write"
                     dt, 0, 0, emptyP, empty, 0, pressure.size(), pressure);
    }

  
  
  MPI_Finalize();
  
  return 0;
}
