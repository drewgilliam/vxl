#ifndef vrel_muse_table_h_
#define vrel_muse_table_h_
//:
//  \file
//  \author Chuck Stewart
//  \date   Summer 2001
//  \date   Modufied May 2004 to store the values in a map
//  \brief  Look-up table for the normalization terms used in the MUSE objective function.
//

#include <iostream>
#include <map>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif

//: Look-up table for the MUSET objective function.
//  Look-up table for the MUSET objective function, derived in James
//  V. Miller's 1997 PhD dissertation at Rensselaer.  An earlier
//  version of appeared in CVPR 1996.  The class computes and stores
//  statistics on the order statistics of Gaussian random variates.
//  Actually, these are for the order statistics of the absolute
//  values of Gaussian random variates.  See vrel_muset_obj for more
//  details.

class vrel_muse_key_type
{
public:
  vrel_muse_key_type(unsigned int k, unsigned int n)
    : k_(k)
    , n_(n)
  {}
  unsigned int k_;
  unsigned int n_;
};

bool
operator<(const vrel_muse_key_type & left_t, const vrel_muse_key_type & right_t);


class vrel_muse_table_entry
{
public:
  vrel_muse_table_entry() = default;
  bool initialized_{ false };
  double expected_;
  double standard_dev_;
  double muse_t_divisor_;
  double muse_t_sq_divisor_;
};

typedef std::map<vrel_muse_key_type, vrel_muse_table_entry> vrel_muse_map_type;


class vrel_muse_table
{
public:
  //: Constructor.
  //  \a table_size is the size of table (= max number of residuals
  //  pre-computed).
  vrel_muse_table(unsigned int /* max_n_stored */) {}

  vrel_muse_table() = default;

  //: Destructor
  ~vrel_muse_table() = default;

  //: Expected value of the kth ordered residual from n samples.
  //  The value is retrieved from the lookup table when possible.
  double
  expected_kth(unsigned int k, unsigned int n);

  //: Standard deviation of the kth ordered residual from n samples.
  //  The value is retrieved from the lookup table when possible.
  double
  standard_dev_kth(unsigned int k, unsigned int n);

  //: The divisor for trimmed statistics.
  //  The value is retrieved from the lookup table when possible.
  double
  muset_divisor(unsigned int k, unsigned int n);


  //: The divisor for trimmed square statistics.
  //  The value is retrieved from the lookup table when possible.
  double
  muset_sq_divisor(unsigned int k, unsigned int n);

private:
  void
  calculate_all(unsigned int k, unsigned int n, vrel_muse_table_entry & entry);

  //: Expected value of the kth ordered residual from n samples.
  //  The value is computed "from scratch".
  double
  calculate_expected(unsigned int k, unsigned int n) const;

  //: Standard deviation of the kth ordered residual from n samples.
  //  The value is computed "from scratch".
  double
  calculate_standard_dev(unsigned int k, unsigned int n, double expected_kth) const;

  //: The divisor for trimmed statistics.
  //  The value is computed "from scratch".
  double
  calculate_divisor(unsigned int k, unsigned int n, double expected_kth) const;

  //: The divisor for trimmed squared statistics.
  //  The value is computed "from scratch".
  double
  calculate_sq_divisor(unsigned int k, unsigned int n, double expected_kth) const;

private:
  vrel_muse_map_type table_;
};

#endif // vrel_muse_table_h_
