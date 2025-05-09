// This is core/vgl/algo/vgl_norm_trans_3d.h
#ifndef vgl_norm_trans_3d_h_
#define vgl_norm_trans_3d_h_
//:
// \file
// \brief The similarity transform that normalizes a 3-d point set
//
// Algorithms to compute projective transformations require that
// the data be conditioned by insuring that the center of gravity
// of the point set is at the origin and the standard deviation
// is isotropic and unity.
//
// \verbatim
//  Modifications
//   Created August 14, 2004 - J.L. Mundy
// \endverbatim

#include <iosfwd>
#include <vnl/vnl_matrix_fixed.h>
#include <vgl/vgl_homg_point_3d.h>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include <vgl/algo/vgl_h_matrix_3d.h>

template <class T>
class vgl_norm_trans_3d : public vgl_h_matrix_3d<T>
{
public:
  // Constructors/Initializers/Destructors-------------------------------------

  vgl_norm_trans_3d();
  vgl_norm_trans_3d(const vgl_norm_trans_3d<T> & M);
  vgl_norm_trans_3d(const vnl_matrix_fixed<T, 4, 4> & M);
  vgl_norm_trans_3d(const T * t_matrix);
  vgl_norm_trans_3d(std::istream & s);
  vgl_norm_trans_3d(const char * filename);
  ~vgl_norm_trans_3d();

  // Operations----------------------------------------------------------------

  //: compute the normalizing transform
  bool
  compute_from_points(const std::vector<vgl_homg_point_3d<T>> & points);

protected: // --- Utility functions -----------------------------------------
  static bool
  scale_xyzroot2(const std::vector<vgl_homg_point_3d<T>> & in, T & radius);

  static void
  center_of_mass(const std::vector<vgl_homg_point_3d<T>> & points, T & cx, T & cy, T & cz);
};

#define VGL_NORM_TRANS_3D_INSTANTIATE(T) extern "please include vgl/algo/vgl_norm_trans_3d.hxx first"

#endif // vgl_norm_trans_3d_h_
