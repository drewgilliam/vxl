#ifndef vrel_homography2d_est_aff_h_
#define vrel_homography2d_est_aff_h_
//:
// \file

#include <vrel/vrel_homography2d_est.h>
#include <vgl/vgl_fwd.h>

//: Subclass of the generalized 8-DOF homography estimator for affine transformations (6 DOF).

class vrel_homography2d_est_aff : public vrel_homography2d_est
{
public:
  //: Constructor from vgl_homg_point_2d's
  vrel_homography2d_est_aff(const std::vector<vgl_homg_point_2d<double>> & from_pts,
                            const std::vector<vgl_homg_point_2d<double>> & to_pts);

  //: Constructor from vnl_vectors
  vrel_homography2d_est_aff(const std::vector<vnl_vector<double>> & from_pts,
                            const std::vector<vnl_vector<double>> & to_pts);

  //: Destructor.
  ~vrel_homography2d_est_aff() override;

  //: Convert a homography to a linear parameter list (for estimation).
  virtual void
  homography_to_parameters(const vnl_matrix<double> & m, vnl_vector<double> & p) const;

  //: Convert a linear parameter list (from estimation) to a homography.
  virtual void
  parameters_to_homography(const vnl_vector<double> & p, vnl_matrix<double> & m) const;
};

#endif // vrel_homography2d_est_aff_h_
