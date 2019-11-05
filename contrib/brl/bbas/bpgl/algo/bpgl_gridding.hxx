// This is bbas/bpgl/algo/bpgl_gridding.hxx
#ifndef bpgl_gridding_hxx_
#define bpgl_gridding_hxx_
//:
// \file
// \brief Transform irregular data to gridded 2D format (e.g. DSMs)
// \author Dan Crispell
// \date Nov 26, 2018
//

#include <iostream>
#include <stdexcept>
#include <memory>
#include <numeric>

#include <vnl/algo/vnl_matrix_inverse.h>
#include <bvgl/bvgl_k_nearest_neighbors_2d.h>
#include "bpgl_gridding.h"


// Inverse distance interpolation
template<class T, class DATA_T>
DATA_T bpgl_gridding::inverse_distance_interp<T, DATA_T>::interp(
    vgl_point_2d<T> interp_loc,
    std::vector<vgl_point_2d<T> > const& neighbor_locs,
    std::vector<DATA_T> const& neighbor_vals,
    T max_dist
    ) const
{
  T weight_sum(0);
  T val_sum(0);
  const unsigned num_neighbors = neighbor_locs.size();
  for (unsigned i=0; i<num_neighbors; ++i) {
    T dist = (neighbor_locs[i] - interp_loc).length();
    if (dist <= max_dist) {
      if (dist < this->dist_eps_) {
        dist = this->dist_eps_;
      }
      T weight = 1.0 / dist;
      weight_sum += weight;
      val_sum += weight*neighbor_vals[i];
    }
  }
  if (weight_sum == T(0)) {
    return this->invalid_val_;
  }
  return val_sum / weight_sum;
}


// Linear interpolation
template<class T, class DATA_T>
DATA_T bpgl_gridding::linear_interp<T, DATA_T>::interp(
    vgl_point_2d<T> interp_loc,
    std::vector<vgl_point_2d<T> > const& neighbor_locs,
    std::vector<DATA_T> const& neighbor_vals,
    T max_dist
    ) const
{
  const unsigned num_neighbors = neighbor_locs.size();

  // vectors of valid neighbor data
  std::vector<double> X,Y,V,W;
  int num_valid_neighbors = 0;

  for (unsigned i=0; i<num_neighbors; ++i) {
    T dist = (neighbor_locs[i] - interp_loc).length();
    if (dist <= max_dist) {
      if (dist < this->dist_eps_) {
        dist = this->dist_eps_;
      }

      // neighbor weight
      double dist_d = static_cast<double>(dist);
      double weight = 1.0 / std::pow(dist_d, dist_iexp_);

      // save to internal storage
      X.emplace_back(static_cast<double>(neighbor_locs[i].x()));
      Y.emplace_back(static_cast<double>(neighbor_locs[i].y()));
      V.emplace_back(static_cast<double>(neighbor_vals[i]));
      W.emplace_back(weight);
      num_valid_neighbors++;
    }
  }

  // check for sufficent neighbors
  if (num_valid_neighbors < 3) {
    // std::cerr << "insufficent neighbors" << std::endl;
    return this->invalid_val_;
  }

  // weight normalization
  double weight_norm = std::accumulate(W.begin(), W.end(), 0.0);
  for (auto& w : W) {
    w /= weight_norm;
  }

  // absolute interpolation: origin at 0
  // relative interpolation: origin at neighbor loc/val centroid
  double x_origin = 0, y_origin = 0, v_origin = 0;
  if (relative_interp_) {
    double N = static_cast<double>(num_valid_neighbors);
    x_origin = std::accumulate(X.begin(), X.end(), 0.0) / N;
    y_origin = std::accumulate(Y.begin(), Y.end(), 0.0) / N;
    v_origin = std::accumulate(V.begin(), V.end(), 0.0) / N;
  }

  // system matrices
  vnl_matrix<double> A(num_valid_neighbors,3);
  vnl_vector<double> b(num_valid_neighbors);
  for (unsigned i=0; i<num_valid_neighbors; ++i) {
    A[i][0] = W[i] * (X[i] - x_origin);
    A[i][1] = W[i] * (Y[i] - y_origin);
    A[i][2] = W[i];
    b[i] = W[i] * (V[i] - v_origin);
  }

  // employ Tikhonov Regularization to cope with degenerate point configurations
  // f = inv(AT * A + lambda * I) * AT * b
  vnl_matrix<double> lambdaI(3, 3, 0);
  for (int d=0; d<3; ++d) {
    lambdaI[d][d] = regularization_lambda_;
  }

  // matrix processing
  vnl_matrix<double> At = A.transpose();
  vnl_matrix<double> AtA_reg = At*A + lambdaI;
  vnl_matrix_inverse<double> inv_AtA_reg(AtA_reg.as_ref());

  // check reciprocal condition number
  auto rcond = inv_AtA_reg.well_condition();
  if (rcond < rcond_thresh_) {
    std::cerr << "matrix has poor condition (" << rcond << ")\n" << std::endl;
    return this->invalid_val_;
  }

  // final solution
  vnl_vector<double> f = inv_AtA_reg.as_matrix() * At * b;
  double x = static_cast<double>(interp_loc.x());
  double y = static_cast<double>(interp_loc.y());
  double value = f[0]*(x - x_origin) + f[1]*(y - y_origin) + f[2] + v_origin;

  // cast as data type
  DATA_T value_return = static_cast<DATA_T>(value);
  return value_return;
}


// interpolate image from scattered locations/values
template<class T, class DATA_T, class INTERP_T>
vil_image_view<DATA_T> bpgl_gridding::grid_data_2d(
    INTERP_T const& interp_fun,
    std::vector<vgl_point_2d<T>> const& data_in_loc,
    std::vector<DATA_T> const& data_in,
    vgl_point_2d<T> out_upper_left,
    size_t out_ni,
    size_t out_nj,
    T step_size,
    unsigned min_neighbors,
    unsigned max_neighbors,
    T max_dist,
    double out_theta_radians
    )
{
  std::cout << "---------------------------------------------------INTERP_FUN TYPE = " << interp_fun.type() << std::endl;

  // total number of points
  size_t npts = data_in_loc.size();

  // validate input
  if (npts != data_in.size()) {
    throw std::runtime_error("Input location and data arrays not equal size");
  }

  // validate min/max neighbor range
  if (size_t(min_neighbors) > npts) {
    throw std::runtime_error("Fewer points than minimum number of neighbors");
  }
  if (size_t(max_neighbors) > npts) {
    max_neighbors = unsigned(npts);
  }
  if (min_neighbors > max_neighbors) {
    throw std::runtime_error("Invalid neighbor range");
  }

  // create knn instance
  bvgl_k_nearest_neighbors_2d<T> knn(data_in_loc);
  if (!knn.is_valid()) {
    throw std::runtime_error("KNN initialization failure");
  }

  vgl_vector_2d<T> i_vec(std::cos(out_theta_radians), std::sin(out_theta_radians));
  vgl_vector_2d<T> j_vec(std::sin(out_theta_radians), -std::cos(out_theta_radians));

  // loop across all grid values
  vil_image_view<DATA_T> gridded(out_ni, out_nj);
  for (unsigned j=0; j<out_nj; ++j) {
    for (unsigned i=0; i<out_ni; ++i) {

      // interpolation point
      vgl_point_2d<T> loc = out_upper_left
                          + i*step_size*i_vec
                          + j*step_size*j_vec;

      // retrieve at most max_neighbors within max_dist of interpolation point
      std::vector<vgl_point_2d<T> > neighbor_locs;
      std::vector<unsigned> neighbor_indices;
      if (!knn.knn(loc, max_neighbors, neighbor_locs, neighbor_indices, max_dist)) {
        throw std::runtime_error("KNN failed to return neighbors");
      }

      // check for at least min_neighbors
      if (neighbor_indices.size() < min_neighbors) {
        gridded(i,j) = interp_fun.invalid_val();
        continue;
      }

      // neighbor values for interpolation
      std::vector<DATA_T> neighbor_vals;
      for (auto nidx : neighbor_indices) {
        neighbor_vals.push_back(data_in[nidx]);
      }

      // interpolate via non-virtual method
      T val = interp_fun(loc, neighbor_locs, neighbor_vals, max_dist);
      gridded(i,j) = val;
    }
  }
  return gridded;
}

// // image gridding
// // template specialization for shared_ptr<base_interp>
// template<class T, class DATA_T>
// vil_image_view<DATA_T>
// grid_data_2d<T, DATA_T, std::shared_ptr<base_interp<T, DATA_T>> >(
//     std::shared_ptr< base_interp<T, DATA_T> > const& interp_ptr,
//     std::vector<vgl_point_2d<T>> const& data_in_loc,
//     std::vector<DATA_T> const& data_in,
//     vgl_point_2d<T> out_upper_left,
//     size_t out_ni,
//     size_t out_nj,
//     T step_size,
//     unsigned min_neighbors = 3,
//     unsigned max_neighbors = 5,
//     T max_dist = vnl_numeric_traits<T>::maxval,
//     double out_theta_radians = 0.0)
// {
//   switch (interp_ptr->type()) {
//     case INVERSE_DISTANCE: {
//       auto inverse_distance_interp_ptr = std::dynamic_pointer_cast< inverse_distance_interp<T, DATA_T> >(interp_ptr);
//       return grid_data_2d(
//           *inverse_distance_interp_ptr, data_in_loc, data_in,
//           out_upper_left, out_ni, out_nj, step_size,
//           min_neighbors, max_neighbors, max_dist, out_theta_radians);
//     }
//     case LINEAR: {
//       auto linear_interp_ptr = std::dynamic_pointer_cast< linear_interp<T, DATA_T> >(interp_ptr);
//       return grid_data_2d(
//           *linear_interp_ptr, data_in_loc, data_in,
//           out_upper_left, out_ni, out_nj, step_size,
//           min_neighbors, max_neighbors, max_dist, out_theta_radians);
//     }
//     default: {
//       throw std::runtime_error("Unrecognized interpolation class")
//     }
//   }
// }


// image to 3D pointset
template<class T, class DATA_T>
vgl_pointset_3d<T> bpgl_gridding::pointset_from_grid(
    vil_image_view<DATA_T> const& grid,
    vgl_point_2d<T> const& upper_left,
    T step_size,
    double out_theta_radians
    )
{
  vgl_pointset_3d<T> ptset;

  vgl_vector_2d<T> i_vec(std::cos(out_theta_radians), std::sin(out_theta_radians));
  vgl_vector_2d<T> j_vec(std::sin(out_theta_radians), -std::cos(out_theta_radians));

  size_t ni = grid.ni(), nj = grid.nj();
  for (size_t j = 0; j<nj; ++j) {
    for (size_t i = 0; i<ni; ++i) {
      if (!std::isfinite(grid(i,j)))
        continue;

      vgl_point_2d<T> pt2d = upper_left + i*step_size*i_vec + j*step_size*j_vec;
      auto pt3d = vgl_point_3d<T>( pt2d.x(), pt2d.y(), static_cast<T>(grid(i,j)) );
      ptset.add_point(pt3d);
    }
  }

  return ptset;
}


// explicit template instantiations macros
#undef BPGL_GRIDDING_INTERP_INSTANIATE
#define BPGL_GRIDDING_INTERP_INSTANIATE(T, DATA_T, INTERP_T) \
template class INTERP_T<T, DATA_T>; \
template vil_image_view<DATA_T> grid_data_2d<T, DATA_T, INTERP_T<T, DATA_T>>( \
    INTERP_T<T, DATA_T> const& interp_fun, \
    std::vector<vgl_point_2d<T>> const& data_in_loc, \
    std::vector<DATA_T> const& data_in, \
    vgl_point_2d<T> out_upper_left, \
    size_t out_ni, \
    size_t out_nj, \
    T step_size, \
    unsigned min_neighbors, \
    unsigned max_neighbors, \
    T max_dist, \
    double out_theta_radians \
    )

#undef BPGL_GRIDDING_INSTANIATE
#define BPGL_GRIDDING_INSTANIATE(T, DATA_T) \
namespace bpgl_gridding { \
BPGL_GRIDDING_INTERP_INSTANIATE(T, DATA_T, inverse_distance_interp); \
BPGL_GRIDDING_INTERP_INSTANIATE(T, DATA_T, linear_interp); \
template vgl_pointset_3d<T> pointset_from_grid<T, DATA_T>( \
    vil_image_view<DATA_T> const& grid, \
    vgl_point_2d<T> const& upper_left, \
    T step_size, \
    double out_theta_radians = 0.0 \
    ); \
}


#endif
