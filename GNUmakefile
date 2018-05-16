HDF5_DIR      	 = /usr/tce/packages/hdf5/hdf5-parallel-1.8.18-intel-18.0.1-openmpi-2.0.0
HDF5INCFLAGS     = -I$(HDF5_DIR)/include
HDF5LIBFLAGS	 = -L$(HDF5_DIR)/lib

MPI_DIR 		 = /usr/tce/packages/openmpi/openmpi-2.0.0-intel-18.0.1
HDFMPIINCFLAGS   = -I$(MPI_DIR)/include
HDFMPILIBFLAGS   = -L$(MPI_DIR)/lib -lhdf5 -lz


CXX = icpc -std=c++11
MPILIBFLAGS=-L/usr/local/lib -lmpi_cxx -lmpi -lm
CXXFLAGS = -g
PPFLAGS = $(HDF5INCFLAGS) $(HDFMPIINCFLAGS)
LDFLAGS = 
LIBFLAGS = $(HDF5LIBFLAGS) $(HDFMPILIBFLAGS) $(MPILIBFLAGS)

all: GEOSDriver.exe ChomboDriver.exe

%.exe: %.o GNUmakefile coupler.o
	$(CXX) $(CXXFLAGS) $(PPFLAGS)  $(LDFLAGS) $< coupler.o $(LIBFLAGS) -o $@

%.o: %.cpp GNUmakefile
	$(CXX) $(CXXFLAGS) -c $(PPFLAGS)  $< -o $@
	$(CXX) -MM $(PPFLAGS) $< > $*.d

vars:
	echo $(CXX)
	echo $(PPFLAGS)
	echo $(LIBFLAGS)
	echo $(LDFLAGS)


clean:
	rm -rf *.hdf5 *.exe *.o *.d *.dSYM

-include $(OBJS:.o=.d)
