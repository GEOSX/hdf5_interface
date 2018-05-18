
#include "coupler.H"

int VCOUNT=9;
int QCOUNT=4;

int quads[] = {0,3,4,1,1,4,5,2,3,6,7,4,4,7,8,5};
double dx=0.5;
double dy=0.5;
double dz=0.0;

inline void createVertex(int globalID, P& p, V& v)
{
  int i = globalID / 3;
  int j = globalID % 3;
  p.X = dx * i;
  p.Y = dy * j;
  p.Z = 0;
  v.X = 0.0001;
  v.Y = -0.0005;
  v.Z = 0.01;
}

//int stop(int, void*){ abort(); return 0;}

int main(int argc, char* argv[])
{

  MPI_Init(&argc, &argv );
  int rank, N;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &N);

  // H5Eset_auto (H5E_DEFAULT,stop, NULL);// for debugging, tell HDF5 to abort on error
  if(rank!= N-1) VCOUNT-=3; //everyone shares 3 vertices, the last rank owns all of it's vertices
  
  // GEOS would replace this logic here.
  int* VID = new int[VCOUNT];   
  P* x = new P[VCOUNT];
  V* v = new V[VCOUNT];
  double* pressure = new double[QCOUNT];
  int* QID = new int[QCOUNT];
 
  for(int i = 0; i < VCOUNT; i++)
  {
    VID[i] = i + 6 * rank; //6 comes from VCOUNT=9 and sharing 3 vertices with neighbor
    createVertex(VID[i], x[i], v[i]);
  }

  for(int i = 0; i < QCOUNT; i++)
  {
    pressure[i] = 0.5;
    QID[i] = i + rank * QCOUNT;
  }
  for(int j = 0; j < QCOUNT * 4; j++)
  {
    quads[j] += 6 * rank;
  }

  int vertexOffset, vertexCount;
  boundaryFileOffsets(MPI_COMM_WORLD, VCOUNT, vertexOffset, vertexCount);
   
  int quadCount, quadOffset;
  boundaryFileOffsets(MPI_COMM_WORLD, QCOUNT, quadOffset, quadCount);

  createBoundaryFile(MPI_COMM_WORLD,
                     "GEOSboundary.hdf5.tmp1", vertexOffset, VCOUNT, vertexCount,
                     VID, x, v, quadOffset, QCOUNT, quadCount, 
                     quads, QID, pressure);

  rename("GEOSboundary.hdf5.tmp1","GEOSboundary.hdf5.tmp2");


  double dt;
  
  rwBoundaryFile(MPI_COMM_WORLD,"GEOSboundary.hdf5.tmp2", true,//true=="read"
                 dt, vertexOffset, VCOUNT, x, v, quadOffset, QCOUNT, pressure);

  V* empty = nullptr;
  P* emptyP = nullptr;
  for(int i = 0; i < QCOUNT; ++i)
  {
    pressure[i] += 0.055;
  }    

  rwBoundaryFile(MPI_COMM_WORLD, "GEOSboundary.hdf5.tmp2", false, //false="write"
                 dt, 0, 0, emptyP, empty, quadOffset, QCOUNT, pressure);

  rename("GEOSboundary.hdf5.tmp2","GEOSboundary.hdf5");
  
  MPI_Finalize();

  delete[] VID;
  delete[] x;
  delete[] v;
  delete[] pressure;
  delete[] QID;
  
  return 0;
}
