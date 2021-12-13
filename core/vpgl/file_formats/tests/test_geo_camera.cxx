#include "testlib/testlib_test.h"

#include <iostream>
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

#include "vpgl/vpgl_lvcs.h"
#include "vpgl/file_formats/vpgl_geo_camera.h"
#include "vpl/vpl.h"
#include "vsl/vsl_binary_io.h"

// helper function: test camera.project(x, y) == (u, v)
void
_test_project(vpgl_geo_camera camera, double x, double y, double u, double v)
{
  double uresult, vresult, uerr, verr;
  camera.project(x, y, 0.0, uresult, vresult);
  std::cout << "project (" << x << "," << y << ") -> "
            << "(" << u << "," << v << "); "
            << "result (" << uresult << "," << vresult << ")\n";
  uerr = std::fabs(u - uresult);
  verr = std::fabs(v - vresult);
  TEST_NEAR("vpgl_geo_camera.project", uerr + verr, 0.0, 1e-3);
}

// helper function test binary io
void
_test_binary_io(vpgl_geo_camera const& cam, std::string filename)
{
  vsl_b_ofstream bp_out(filename.c_str());
  TEST("Created " + filename + "for writing", (!bp_out), false);

  cam.b_write(bp_out);
  bp_out.close();

  vsl_b_ifstream bp_in(filename.c_str());
  TEST("Opened " + filename + " for reading", (!bp_in), false);

  vpgl_geo_camera cam_read;
  cam_read.b_read(bp_in);
  bp_in.close();

  std::cout << "Recovered vpgl_geo_camera:\n" << cam_read << std::endl;
  TEST("Recovery from binary read", cam_read, cam);

  // remove file:
  vpl_unlink(filename.ctr());
}


static void
test_geo_camera()
{
  double x, y, u, v, lx, ly, lz;
  std::vector<std::tuple<double, double, double, double> > data;

  // WGS84 geotransform
  // https://gdal.org/user/raster_data_model.html#affine-geotransform
  // Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  // Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  std::array<double, 6> geotransform = {100.0, 1e-6, 0.0, 30.0, 0.0, -1e-6};
  std::cout << "geotransform for testing: ";
  for (auto const& v : geotransform)
    std::cout << v << " ";
  std::cout << std::endl;

  // WGS84 geocam with no lvcs
  // global coordinate -> raster coordinate
  auto cam0 = load_geo_camera_from_geotransform(geotransform);
  std::cout << "geo_camera: WGS84 without LVCS:\n" << cam0 << std::endl;

  data.clear();
  data.push_back(std::make_tuple(100.0, 30.0, 0.0, 0.0));
  data.push_back(std::make_tuple(100.0 + (50*1e-6), 30.0 - (50*1e-6), 50.0, 50.0));

  for (auto const& item : data) {
    std::tie(x, y, u, v) = item;
    _test_project(cam0, x, y, u, v);
  }

  _test_binary_io(cam0, "test_geo_camera_io.tmp");

  // WGS84 geocam with lvcs
  // local coordinate -> raster coordinate
  auto lvcs = vpgl_lvcs(30.0, 100.0, 0.0);
  auto cam1 = load_geo_camera_from_geotransform(geotransform, -1, 0, &lvcs);
  std::cout << "geo_camera: WGS84 with LVCS:\n" << cam1 << "\n";

  data.clear();
  data.push_back(std::make_tuple(0.0, 0.0, 0.0, 0.0));

  lvcs.global_to_local(100.0 + (50*1e-6), 30.0 - (50*1e-6), 0.0,
                       vpgl_lvcs::wgs84, lx, ly, lz);
  data.push_back(std::make_tuple(lx, ly, 50.0, 50.0));

  for (auto const& item : data) {
    std::tie(x, y, u, v) = item;
    _test_project(cam1, x, y, u, v);
  }

  _test_binary_io(cam1, "test_geo_camera_lvcs_io.tmp");

}

TESTMAIN(test_geo_camera);
