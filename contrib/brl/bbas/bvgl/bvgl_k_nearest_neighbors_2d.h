// This is brl/bbas/bvgl/bvgl_k_nearest_neighbors_2d.h
#ifndef bvgl_k_nearest_neighbors_2d_h_
#define bvgl_k_nearest_neighbors_2d_h_
//:
// \file
// \brief Uses the nabo knn algorithm to find nearest neighbors
// \author Daniel Crispell
//

#include <iostream>
#include <iosfwd>
#include <limits>
#include <algorithm>
#include <utility>
#include <memory>
#include <stdexcept>
#include <bnabo/bnabo.h>
#include <vgl/vgl_point_2d.h>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif

template <class Type>
class bvgl_k_nearest_neighbors_2d
{
public:
  //: Construct from set of points
  bvgl_k_nearest_neighbors_2d(std::vector<vgl_point_2d<Type> > const &ptset, Type tolerance = Type(0));

  //: check class validity
  bool is_valid() const;

  //: query for the closest point (k = 1) including self
  bool closest_point(vgl_point_2d<Type> const& p, vgl_point_2d<Type>& cp) const;

  //: query for the index in the pointset with the closest point including self
  bool closest_index(vgl_point_2d<Type> const& p, unsigned& index) const;

  //: find the indices of the k closest neighbors (limited by max_dist)
  bool knn_indices(
      vgl_point_2d<Type> const& p,
      unsigned k,
      std::vector<unsigned> &neighbor_indices,
      Type max_dist = std::numeric_limits<Type>::infinity()
      ) const;

  //: find k nearest neighbors
  bool knn(
      vgl_point_2d<Type> const& p,
      unsigned k,
      std::vector<vgl_point_2d<Type>>& neighbor_locs,
      Type max_dist = std::numeric_limits<Type>::infinity()
      ) const;

  //: find k nearest neighbors and indices
  bool knn(
      vgl_point_2d<Type> const& p,
      unsigned k,
      std::vector<vgl_point_2d<Type>>& neighbor_locs,
      std::vector<unsigned>& neighbor_indices,
      Type max_dist = std::numeric_limits<Type>::infinity()
      ) const;

protected:
  //: creates and populates the underlying data structure
  bool create();

  Type tolerance_;
  std::unique_ptr<Nabo::NearestNeighbourSearch<Type>> search_tree_;
  vnl_matrix<Type> M_;//a matrix form(2 x n) of the pointset used by nabo
  std::vector<vgl_point_2d<Type>> ptset_;
  unsigned flags_;//control various actions during queries
};


template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::create(){
  flags_ = 0;
  flags_ = flags_ |  Nabo::NearestNeighbourSearch<Type>::ALLOW_SELF_MATCH;
  unsigned n = ptset_.size(), dim = 2;
  if(n==0){
    search_tree_ = nullptr;
    return false;
  }
  M_.set_size(dim, n);
  for(unsigned i = 0; i<n; ++i){
    vgl_point_2d<Type> pi = ptset_[i];
    M_[0][i]=pi.x();
    M_[1][i]=pi.y();
  }
  search_tree_.reset(Nabo::NearestNeighbourSearch<Type>::createKDTreeLinearHeap(M_, dim));
  return true;
}


template <class Type>
bvgl_k_nearest_neighbors_2d<Type>::bvgl_k_nearest_neighbors_2d(
    std::vector<vgl_point_2d<Type>> const& ptset,
    Type tolerance):
  tolerance_(tolerance),
  ptset_(ptset)
{
  create();
}

template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::is_valid() const
{
  if(!search_tree_)
    return false;
  else
    return true;
}

// main workhorse function:
// discover k-nearest neighbors to point p within max_dist radius
template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::knn_indices(
    vgl_point_2d<Type> const& p,
    unsigned k,
    std::vector<unsigned> &neighbor_indices,
    Type max_dist
    ) const
{
  // clear output
  neighbor_indices.clear();

  // check search tree validity
  if(!search_tree_) {
    return false;
  }

  // variable init
  vnl_vector<Type> q(2), dists2(k);
  vnl_vector<int> indices(k);
  q[0]=p.x();
  q[1]=p.y();

  // search
  search_tree_->knn(q, indices, dists2, k, tolerance_, flags_, max_dist);

  // gather output
  for (unsigned i = 0; i < k; ++i) {

    // infinite distance or invalid index == no more neighbors found
    // if max_dist is finite, fewer than k neighbors is fine
    // if max_dist is infinite, fewer than k neighbors is unexpected
    if(!std::isfinite(dists2[i]) || indices[i]<0) {
      if (std::isfinite(max_dist)) {
        break;
      } else {
        return false;
      }
    }

    // add valid index to list
    neighbor_indices.push_back(static_cast<unsigned>(indices[i]));
  }
}


template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::closest_index(
    vgl_point_2d<Type> const& p,
    unsigned& ci
    ) const
{
  std::vector<unsigned> indices;
  if (!this->knn_indices(p, 1, indices))
    return false;
  ci = indices[0];
  return true;
}


template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::closest_point(
    vgl_point_2d<Type> const& p,
    vgl_point_2d<Type>& cp
    ) const
{
  unsigned index = 0;
  if (!this->closest_index(p, index))
    return false;
  cp = ptset_[index];
  return true;
}


template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::knn(
    vgl_point_2d<Type> const& p,
    unsigned k,
    std::vector<vgl_point_2d<Type>>& neighbor_locs,
    Type max_dist
    ) const
{
  std::vector<unsigned> neighbor_indices;
  return this->knn(p, k, neighbor_locs, neighbor_indices, max_dist);
}


template <class Type>
bool bvgl_k_nearest_neighbors_2d<Type>::knn(
    vgl_point_2d<Type> const& p,
    unsigned k,
    std::vector<vgl_point_2d<Type>>& neighbor_locs,
    std::vector<unsigned>& neighbor_indices,
    Type max_dist
    ) const
{
  neighbor_locs.clear();
  if (!this->knn_indices(p, k, neighbor_indices, max_dist))
    return false;
  for(auto i : neighbor_indices)
    neighbor_locs.push_back(ptset_[i]);
  return true;
}

#endif // bvgl_k_nearest_neighbors_2d_h_
