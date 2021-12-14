#ifndef vpgl_io_geo_camera_h_
#define vpgl_io_geo_camera_h_
//:
// \file
#include <vsl/vsl_binary_io.h>
#include <vpgl/file_formats/vpgl_geo_camera.h>

//: Binary save lvcs to stream
void vsl_b_write(vsl_b_ostream & os, vpgl_geo_camera const& cam);

//: Binary load lvcs from stream.
void vsl_b_read(vsl_b_istream & is, vpgl_geo_camera &cam);

//: Print human readable summary of object to a stream
void vsl_print_summary(std::ostream& os, vpgl_geo_camera const& cam);

#endif
