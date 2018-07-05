# Define a variable basic_geometry_source_files containing
# a list of the source files for the basic_geometry library
set(basic_geometry_source_files
    PorousWall.cpp
)

# Create a library called basic_boundary_conditions which includes the 
# source files defined in basic_boundary_conditions_source_files
add_library(basic_geometry ${basic_geometry_source_files})

TARGET_LINK_LIBRARIES(basic_geometry
  SAMRAI_appu SAMRAI_geom SAMRAI_solv SAMRAI_algs SAMRAI_mesh
  SAMRAI_math SAMRAI_pdat SAMRAI_xfer SAMRAI_hier SAMRAI_tbox
  ${HDF5_C_LIBRARIES} gfortran)

if (WITH_ZLIB)
    TARGET_LINK_LIBRARIES(basic_geometry z)
endif()
