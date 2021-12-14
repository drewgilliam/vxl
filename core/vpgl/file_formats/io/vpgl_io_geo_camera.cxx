#include "vpgl/file_formats/io/vpgl_io_geo_camera.h"
#include "vpgl/io/vpgl_io_lvcs.h"
#include "vnl/io/vnl_io_matrix.h"

//: Binary write to stream.
void
vsl_b_write(vsl_b_ostream & os, vpgl_geo_camera const& cam)
{
  if (!os)
    return;
  unsigned version = 2;
  vsl_b_write(os, version);

  // transformation matrix
  vsl_b_write(os, cam.trans_matrix());

  // write LVCS i available
  bool has_lvcs = cam.has_lvcs();
  vsl_b_write(os, has_lvcs);
  if (has_lvcs) {
    vsl_b_write(os, cam.lvcs());
  }

  // additional info
  vsl_b_write(os, cam.scale_tag());
  vsl_b_write(os, cam.utm_zone());
  vsl_b_write(os, cam.south_flag());

  double sx, sy;
  cam.pixel_spacing(sx, sy);
  vsl_b_write(os, sx);
  vsl_b_write(os, sy);

}

//: Binary read from stream.
void
vsl_b_read(vsl_b_istream & is, vpgl_geo_camera & cam)
{
  if (!is)
    return;
  short ver;
  vsl_b_read(is, ver);
  switch (ver)
  {
    case 1: {
      unsigned nrows, ncols;
      vsl_b_read(is, nrows);
      vsl_b_read(is, ncols);

      vnl_matrix<double> trans_matrix(nrows, ncols, 0.0);
      for (unsigned i = 0; i < nrows; i++)
        for (unsigned j = 0; j < ncols; j++)
          vsl_b_read(is, trans_matrix[i][j]);

      vpgl_lvcs_sptr lvcs;
      vsl_b_read(is, lvcs);

      bool is_utm;
      vsl_b_read(is, is_utm);
      int utm_zone;
      vsl_b_read(is, utm_zone);
      int south_flag;
      vsl_b_read(is, south_flag);
      bool scale_tag;
      vsl_b_read(is, scale_tag);

      vpgl_geo_camera temp(trans_matrix, lvcs, scale_tag,
                           utm_zone, south_flag);
      cam = temp;
      break;
    }
    case 2: {
      vnl_matrix<double> trans_matrix;
      vsl_b_read(is, trans_matrix);

      bool has_lvcs;
      vpgl_lvcs_sptr lvcs = nullptr;
      vsl_b_read(is, has_lvcs);
      if (has_lvcs) {
        vsl_b_read(is, lvcs);
      }

      bool scale_tag;
      vsl_b_read(is, scale_tag);
      int utm_zone;
      vsl_b_read(is, utm_zone);
      bool south_flag;
      vsl_b_read(is, south_flag);

      double sx, sy;
      vsl_b_read(is, sx);
      vsl_b_read(is, sy);

      vpgl_geo_camera temp(trans_matrix, lvcs, scale_tag,
                           utm_zone, south_flag, sx, sy);
      cam = temp;
      break;
    }
    default:
      std::cerr << "I/O ERROR: vpgl_geo_camera::b_read(vsl_b_istream&)\n"
                << "           Unknown version number " << ver << '\n';
      is.is().clear(std::ios::badbit); // Set an unrecoverable IO error on stream
      return;
  }
}

//: Print human readable summary of object to a stream
void
vsl_print_summary(std::ostream & os, vpgl_geo_camera const& cam)
{
  os << cam << '\n';
}
