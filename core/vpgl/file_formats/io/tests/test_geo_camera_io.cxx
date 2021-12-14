#include <iostream>
#include <cmath>
#include "testlib/testlib_test.h"
#ifdef _MSC_VER
#  include "vcl_msvc_warnings.h"
#endif
#include "vpgl/file_formats/vpgl_geo_camera.h"
#include "vpgl/file_formats/io/vpgl_io_geo_camera.h"
#include "vpl/vpl.h"
#include "vsl/vsl_binary_io.h"

void
_test_binary_io(vpgl_geo_camera const& cam, std::string filename,
                std::string title)
{

  std::cout << "\n===== " << title << " =====\n";
  std::cout << "using " << filename << "\n";
  std::cout << "\ntest vpgl_geo_camera:\n" << cam << "\n";

  vsl_b_ofstream os(filename.c_str());
  TEST("Created file for writing", (!os), false);

  vsl_b_write(os, cam);
  os.close();
  TEST("Wrote vpgl_geo_camera to file", true, true);

  vsl_b_ifstream is(filename.c_str());
  TEST("Opened file for reading", (!is), false);

  vpgl_geo_camera cam_read;
  vsl_b_read(is, cam_read);
  is.close();

  std::cout << "\nRecovered vpgl_geo_camera:\n" << cam_read << std::endl;
  TEST("Recovery from binary read", cam_read, cam);

  vpl_unlink(filename.c_str());
}

void
_test_wgs84()
{
  // WGS84 geotransform
  // https://gdal.org/user/raster_data_model.html#affine-geotransform
  // Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  // Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  std::array<double, 6> geotransform = {100.0, 1e-6, 0.0, 30.0, 0.0, -1e-6};
  std::cout << "WGS84 geotransform for testing: ";
  for (auto const& v : geotransform)
    std::cout << v << " ";
  std::cout << std::endl;

  // WGS84 lvcs
  vpgl_lvcs lvcs(30.0, 100.0, 0.0);

  // WGS84 geocam with no lvcs
  // global coordinate -> raster coordinate
  auto cam0 = load_geo_camera_from_geotransform(geotransform);
  _test_binary_io(cam0, "test_geo_camera_io_wgs84.tmp",
                  "WGS84 without LVCS");

  // WGS84 geocam with lvcs
  // local coordinate -> raster coordinate
  auto cam1 = load_geo_camera_from_geotransform(geotransform, -1, 0, &lvcs);
  _test_binary_io(cam1, "test_geo_camera_io_wgs84_lvcs.tmp",
                  "WGS84 with LVCS");
}

void
_test_utm()
{
  // UTM location
  double lat = 38.982859, lon = -117.057278, elev = 1670;
  double easting = 495039.001, northing = 4314875.991;
  int utm_zone = 11;
  bool south_flag = 0;

  // UTM geotransform
  std::array<double, 6> geotransform = {easting, 0.5, 0.0, northing, 0.0, -0.5};
  std::cout << "UTM geotransform for testing: ";
  for (auto const& v : geotransform)
    std::cout << v << " ";
  std::cout << std::endl;

  // UTM lvcs
  vpgl_lvcs lvcs(lat, lon, elev, vpgl_lvcs::utm,
                 vpgl_lvcs::DEG, vpgl_lvcs::METERS);

  // UTM geocam with no lvcs
  // global coordinate -> raster coordinate
  auto cam0 = load_geo_camera_from_geotransform(geotransform, utm_zone, south_flag);
  _test_binary_io(cam0, "test_geo_camera_io_utm.tmp",
                  "UTM without LVCS");

  // UTM geocam with lvcs
  // local coordinate -> raster coordinate
  auto cam1 = load_geo_camera_from_geotransform(geotransform, utm_zone, south_flag, &lvcs);
  _test_binary_io(cam1, "test_geo_camera_io_utm_lvcs.tmp",
                  "UTM with LVCS");

}

static void
test_geo_camera_io()
{
  _test_wgs84();
  _test_utm();
}

TESTMAIN(test_geo_camera_io);
