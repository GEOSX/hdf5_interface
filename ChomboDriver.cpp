#include "coupler.H"

int stop(int, void*){ abort(); return 0;}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv );
  int rank, N;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &N);

  H5Eset_auto (H5E_DEFAULT,stop, NULL);

#if 0

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

  dt = 1.155;
  if (rank==0)
  {
    writeBoundaryFile(MPI_COMM_NULL, "GEOSboundary.hdf5", dt, 0, VCOUNT, x, v,
                      0, QCOUNT, pressure);
  }
#endif

  MPI_Finalize();
  return 0;
}
