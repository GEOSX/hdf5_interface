HDF5_DIR      = /Users/bvs/hdf5-1.8.18
HDFMPIINCFLAGS   = -I$(HDF5_DIR).parallel/include
HDFMPILIBFLAGS   = -L$(HDF5_DIR).parallel/lib -lhdf5  -lz


CXX = clang++ -std=c++11
MPILIBFLAGS=-L/usr/local/lib -lmpi_cxx -lmpi -lm
CXXFLAGS = -g
PPFLAGS = $(HDFMPIINCFLAGS)
LDFLAGS = 
LIBFLAGS = $(HDFMPILIBFLAGS) $(MPILIBFLAGS)

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
