# core/vpgl/file_formats/io/Cmakelists.txt

set(vpgl_file_formats_io_sources
    vpgl_io_geo_camera.h              vpgl_io_geo_camera.cxx
)

vxl_add_library(LIBRARY_NAME ${VXL_LIB_PREFIX}vpgl_file_formats_io LIBRARY_SOURCES ${vpgl_file_formats_io_sources})

target_link_libraries(${VXL_LIB_PREFIX}vpgl_file_formats_io ${VXL_LIB_PREFIX}vpgl_file_formats ${VXL_LIB_PREFIX}vpgl_io ${VXL_LIB_PREFIX}vsl)

if( BUILD_TESTING )
  add_subdirectory(tests)
endif()
