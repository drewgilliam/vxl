#include "vpgl_io_lvcs.h"
//:
// \file
#include "vpgl/vpgl_lvcs.h"
#include <vnl/io/vnl_io_matrix_fixed.h>

void
vsl_b_write(vsl_b_ostream & os, vpgl_lvcs const & lvcs)
{
  if (!os)
  {
    return;
  }

  unsigned version = 2;
  vsl_b_write(os, version);

  auto csn = static_cast<unsigned>(lvcs.get_cs_name());
  vsl_b_write(os, csn);

  double lat, lon, elev;
  lvcs.get_origin(lat, lon, elev);
  vsl_b_write(os, lat);
  vsl_b_write(os, lon);
  vsl_b_write(os, elev);

  double lat_scale, lon_scale;
  lvcs.get_scale(lat_scale, lon_scale);
  vsl_b_write(os, lat_scale);
  vsl_b_write(os, lon_scale);

  auto gaunit = static_cast<unsigned>(lvcs.geo_angle_unit());
  vsl_b_write(os, gaunit);

  auto xyzunit = static_cast<unsigned>(lvcs.local_length_unit());
  vsl_b_write(os, xyzunit);

  double lox, loy, theta;
  lvcs.get_transform(lox, loy, theta);
  vsl_b_write(os, lox);
  vsl_b_write(os, loy);
  vsl_b_write(os, theta);

  double easting, northing, utm_elev;
  int utm_zone;
  bool south_flag;
  lvcs.get_utm_origin(easting, northing, utm_elev, utm_zone, south_flag);
  vsl_b_write(os, utm_zone);
  vsl_b_write(os, south_flag);

}

//: Binary load lvcs from stream.
void
vsl_b_read(vsl_b_istream & is, vpgl_lvcs & lvcs)
{
  if (!is)
  {
    return;
  }

  short ver;
  vsl_b_read(is, ver);

  switch (ver)
  {
    case 1: {
      unsigned cs_name;
      vsl_b_read(is, cs_name);

      double lat, lon, elev;
      vsl_b_read(is, lat);
      vsl_b_read(is, lon);
      vsl_b_read(is, elev);

      double lat_scale, lon_scale;
      vsl_b_read(is, lat_scale);
      vsl_b_read(is, lon_scale);

      unsigned ang_unit, len_unit;
      vsl_b_read(is, ang_unit);
      vsl_b_read(is, len_unit);

      double lox, loy, theta;
      vsl_b_read(is, lox);
      vsl_b_read(is, loy);
      vsl_b_read(is, theta);

      vpgl_lvcs temp(lat, lon, elev,
                     static_cast<vpgl_lvcs::cs_names>(cs_name),
                     lat_scale, lon_scale,
                     static_cast<vpgl_lvcs::AngUnits>(ang_unit),
                     static_cast<vpgl_lvcs::LenUnits>(len_unit),
                     lox, loy, theta);
      lvcs = temp;
      break;
    }
    case 2: {
      unsigned cs_name;
      vsl_b_read(is, cs_name);

      double lat, lon, elev;
      vsl_b_read(is, lat);
      vsl_b_read(is, lon);
      vsl_b_read(is, elev);

      double lat_scale, lon_scale;
      vsl_b_read(is, lat_scale);
      vsl_b_read(is, lon_scale);

      unsigned ang_unit, len_unit;
      vsl_b_read(is, ang_unit);
      vsl_b_read(is, len_unit);

      double lox, loy, theta;
      vsl_b_read(is, lox);
      vsl_b_read(is, loy);
      vsl_b_read(is, theta);

      int utm_zone;
      bool south_flag;
      vsl_b_read(is, utm_zone);
      vsl_b_read(is, south_flag);

      vpgl_lvcs temp(lat, lon, elev,
                     static_cast<vpgl_lvcs::cs_names>(cs_name),
                     lat_scale, lon_scale,
                     static_cast<vpgl_lvcs::AngUnits>(ang_unit),
                     static_cast<vpgl_lvcs::LenUnits>(len_unit),
                     lox, loy, theta,
                     utm_zone, south_flag);
      lvcs = temp;
      break;
    }
    default:
      std::cerr << "I/O ERROR: vpgl_lvcs::b_read(vsl_b_istream&)\n"
                << "           Unknown version number " << ver << '\n';
      is.is().clear(std::ios::badbit); // Set an unrecoverable IO error on stream
      return;
  }
}

//: Print human readable summary of object to a stream
void
vsl_print_summary(std::ostream & os, const vpgl_lvcs & c)
{
  os << c << '\n';
}

//: Binary save lvcs sptr to stream
void
vsl_b_write(vsl_b_ostream & os, vpgl_lvcs_sptr const & lvcs_sptr)
{
  if (!lvcs_sptr)
    return;
  vpgl_lvcs * lvcs = lvcs_sptr.ptr();
  vsl_b_write(os, *lvcs);
}

//: Binary load lvcs sptr from stream.
void
vsl_b_read(vsl_b_istream & is, vpgl_lvcs_sptr & lvcs_sptr)
{
  vpgl_lvcs * lvcs = new vpgl_lvcs();
  vsl_b_read(is, *lvcs);
  lvcs_sptr = lvcs;
}
