// This is bbas/bpgl/algo/bpgl_gridding.h
#ifndef bpgl_gridding_h_
#define bpgl_gridding_h_
//:
// \file
// \brief Transform irregular data to gridded 2D format (e.g. DSMs)
// \author Dan Crispell
// \date Nov 26, 2018
//

#include <vector>
#include <limits>

#include <vil/vil_image_view.h>
#include <vgl/vgl_pointset_3d.h>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif

namespace bpgl_gridding
{


//: Interpolation abstract class
// invalid_val: Value to return when interpolation is not appropriate
// dist_eps: The smallest meaningful distance of input points. Must be > 0.
template<class T, class DATA_T>
class base_interp
{
 public:

  // constructors
  base_interp() = default;

  base_interp(
      DATA_T invalid_val,
      T dist_eps) :
    invalid_val_(invalid_val),
    dist_eps_(dist_eps)
  {}

  // accessors
  DATA_T invalid_val() const { return invalid_val_; }
  void invalid_val(DATA_T x) { invalid_val_ = x; }

  T dist_eps() const { return dist_eps_; }
  void dist_eps(T x) { dist_eps_ = x; }

  // interpolation operator
  virtual DATA_T operator() (
      vgl_point_2d<T> interp_loc,
      std::vector<vgl_point_2d<T> > const& neighbor_locs,
      std::vector<DATA_T> const& neighbor_vals,
      T max_dist
      ) const = 0;

 protected:

  // parameters with defaults
  DATA_T invalid_val_ = DATA_T(NAN);
  T dist_eps_ = 1e-5;

};


//: Inverse distance interpolation class
// base class parameters (see "base_interp" above)
template<class T, class DATA_T>
class inverse_distance_interp : public base_interp<T, DATA_T>
{
 public:

  // constructors (inherit from base_interp)
  using base_interp<T, DATA_T>::base_interp;

  // virtual interpolation operator
  DATA_T operator() (
      vgl_point_2d<T> interp_loc,
      std::vector<vgl_point_2d<T> > const& neighbor_locs,
      std::vector<DATA_T> const& neighbor_vals,
      T max_dist = std::numeric_limits<T>::infinity()
      ) const override;

};


//: Linear interpolation class
// base class parameters (see "base_interp" above)
// dist_iexp: neighbor weight is proportional to (1/dist)^dist_iexp
// regularization_lambda: Larger regularization values will bias the solution
//    towards "flatter" functions.  Very large values will result in weighted
//    averages of neighbor values.
// rcond_thresh: Threshold for inverse matrix conditioning. Must be > 0
// relative_interp: interpolation relative to neighbor centroids
//
// Note: internals are represented at double precision to ensure
// accuracy of the final result.  Only on output is the resulting
// interpolated value cast back to DATA_T
template<class T, class DATA_T>
class linear_interp : public base_interp<T, DATA_T>
{
 public:

  // constructors (inherit from base_interp)
  using base_interp<T, DATA_T>::base_interp;

  // accessors
  int dist_iexp() const { return dist_iexp_; }
  void dist_iexp(int x) { dist_iexp_ = x; }

  double regularization_lambda() const { return regularization_lambda_; }
  void regularization_lambda(double x) { regularization_lambda_ = x; }

  double rcond_thresh() const { return rcond_thresh_; }
  void rcond_thresh(double x) { rcond_thresh_ = x; }

  bool relative_interp() const { return relative_interp_; }
  void relative_interp(bool x) { relative_interp_ = x; }

  // interpolation operator
  DATA_T operator() (
      vgl_point_2d<T> interp_loc,
      std::vector<vgl_point_2d<T> > const& neighbor_locs,
      std::vector<DATA_T> const& neighbor_vals,
      T max_dist = std::numeric_limits<T>::infinity()
      ) const override;

 private:

  // parameters with defaults
  int dist_iexp_ = 2;
  double regularization_lambda_ = 1e-3;
  double rcond_thresh_ = 1e-8;
  bool relative_interp_ = true;

};


//: image interpolation from scattered locations/values
template<class T, class DATA_T>
vil_image_view<DATA_T> grid_data_2d(
    base_interp<T, DATA_T> const& interp_fun,
    std::vector<vgl_point_2d<T>> const& data_in_loc,
    std::vector<DATA_T> const& data_in,
    vgl_point_2d<T> out_upper_left,
    size_t out_ni,
    size_t out_nj,
    T step_size,
    unsigned min_neighbors = 3,
    unsigned max_neighbors = 5,
    T max_dist = std::numeric_limits<T>::infinity(),
    double out_theta_radians = 0.0
    );

// image as 3D point cloud
template<class T, class DATA_T>
vgl_pointset_3d<T> pointset_from_grid(
    vil_image_view<DATA_T> const& grid,
    vgl_point_2d<T> const& upper_left,
    T step_size,
    double out_theta_radians = 0.0
    );


} // end namespace bpgl_gridding

#endif
