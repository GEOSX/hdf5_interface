#
# Specify all headers
#
set(hdf5_interface_headers
  coupler.hpp )

#
# Specify all sources
#
set(hdf5_interface_sources
  coupler.cpp )

set( dependencyList hdf5 mpi )

blt_add_library( NAME                  hdf5_interface
                 SOURCES               ${hdf5_interface_sources}
                 HEADERS               ${hdf5_interface_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                TRUE
               )
               
target_include_directories( hdf5_interface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX hdf5_interface )
