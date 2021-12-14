add_executable( vpgl_file_formats_io_test_all
  test_driver.cxx
  test_geo_camera_io.cxx
)

target_link_libraries( vpgl_file_formats_io_test_all ${VXL_LIB_PREFIX}vpgl_file_formats_io ${VXL_LIB_PREFIX}vpl ${VXL_LIB_PREFIX}testlib )

add_test( NAME vpgl_test_geo_camera_io COMMAND $<TARGET_FILE:vpgl_file_formats_io_test_all> test_geo_camera_io )

add_executable( vpgl_file_formats_io_test_include test_include.cxx )
target_link_libraries( vpgl_file_formats_io_test_include ${VXL_LIB_PREFIX}vpgl_file_formats_io )
