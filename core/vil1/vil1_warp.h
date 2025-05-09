// This is core/vil1/vil1_warp.h
#ifndef vil1_warp_h_
#define vil1_warp_h_
//:
// \file
// \author awf@robots.ox.ac.uk
// \date 04 Dec 00

#include "vil1_image.h"
#include "vil1_memory_image_of.h"

//: Class to encapsulate 2d-2d mapping.
class vil1_warp_mapping
{
public:
  virtual ~vil1_warp_mapping() = default;
  virtual void
  forward_map(double x1, double y1, double * x2, double * y2) const = 0;
  virtual void
  inverse_map(double x2, double y2, double * x1, double * y1) const = 0;
};

//: Enum which selects type of interpolation for vil1_warp*
enum vil1_warp_interpolation_type
{
  vil1_warp_interpolation_nearest_neighbour,
  vil1_warp_interpolation_bilinear,
  vil1_warp_interpolation_bicubic
};

//: Warp an image under a 2D map.
// The output image is set to the same size as the input image if the
// output width and height arguments are omitted or set to negative
// values.
// The mapping is such that out(mapper.forward_map(x,y)) = in(x,y);
vil1_image
vil1_warp(const vil1_image & in, const vil1_warp_mapping & mapper, int out_width = -1, int out_height = -1);

//: Warp, specifying interpolation
vil1_image
vil1_warp(const vil1_image & in,
          const vil1_warp_mapping & mapper,
          vil1_warp_interpolation_type i,
          int out_width = -1,
          int out_height = -1);

//: Templated warper
template <class PixelType, class Mapper>
void
vil1_warp_output_driven(const vil1_memory_image_of<PixelType> & in,
                        vil1_memory_image_of<PixelType> & out,
                        const Mapper & mapper,
                        vil1_warp_interpolation_type);

#endif // vil1_warp_h_
