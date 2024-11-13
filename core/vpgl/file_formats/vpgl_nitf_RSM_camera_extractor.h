// This is core/vpgl/file_formats/vpgl_nitf_RSM_camera_extractor.h
#ifndef vpgl_nitf_RSM_camera_extractor_h_
#define vpgl_nitf_RSM_camera_extractor_h_
//:
// \file
// \brief: instance a replacement sensor model (RSM) polynomial camera
// \author J.L. Mundy
// \date October 2023
//
#include <string>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include <vpgl/vpgl_RSM_camera.h>
#include <vil/file_formats/vil_nitf2_image.h>
#include <vnl/vnl_double_2.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_box_3d.h>
#include <iostream>
#include <vgl/vgl_polygon.h>
#include <vpgl_replacement_sensor_model_tres.h>
#include <vnl/vnl_matrix_fixed.h>
struct image_time {
  int year, month, day, hour, min, sec;
};
// information extracted from NITF header
struct rsm_metadata{
  std::string catalog_id_;                     bool catalog_id_valid=false;
  std::string platform_name_;                  bool platform_name_valid=false;
  std::string image_name_;                     bool image_name_valid = false;
  std::vector<int> acquisition_time_;          bool acquisition_time_valid = false;
  unsigned effective_bits_per_pixel_;          bool effective_bits_per_pixel_valid = false;
  std::string image_type_;                     bool image_type_valid = false;
  std::string igeolo_;                         bool igeolo_valid = false;
  bool xy_corners_valid = false;
  bool xyz_corners_valid = false;
  vgl_point_2d<double> upper_left_;            
  vgl_point_2d<double> upper_right_;
  vgl_point_2d<double> lower_left_;
  vgl_point_2d<double> lower_right_;
  vgl_box_3d<double> bounding_box_; 
  vgl_polygon<double> footprint_;
  vgl_point_2d<double> image_offset_;           bool image_offset_valid = false;
  vgl_point_2d<double> rsm_image_offset_;       bool rsm_image_offset_valid = false;
  //       x => column y => row
  vgl_point_2d<double> min_image_corner_;       bool image_corners_valid = false;
  vgl_point_2d<double> max_image_corner_;
  bool geodetic_=true;
  //vs. 0->360
  bool longitude_plus_minus_180_=true; 
  bool rectangular_=false;
  vgl_point_3d<double> rect_origin_;            bool rect_origin_valid = false;
  vnl_matrix_fixed<double, 3, 3> rect_trans_;   bool rect_trans_valid = false;
  double sun_azimuth_radians_=0.0;              bool sun_azimuth_valid = false;
  double sun_elevation_radians_=0.0;            bool sun_elevation_valid = false;
};
// if image is cropped defines offset
struct ichipb_data{
  bool ichipb_data_valid = false;
  std::pair<double, double> translation_;
  std::vector<std::pair<double, double> > F_grid_points_;
  std::vector<std::pair<double, double> > O_grid_points_;
  double scale_factor_;
  bool anamorphic_corr_;
};
struct adjustable_parameter_covar{
  bool defined_ = false;
  size_t num_adj_params_ = 0;
  size_t num_original_adj_params_ = 0;

  // local coordinate system
  vnl_vector_fixed<double, 3> translation_;
  vnl_matrix_fixed<double, 3, 3> rotation_;
  
  //  adjustable param    index  (1-based)
  std::map<std::string, int> covar_index_;

  //   group id         covariance matrix
  std::vector<vnl_matrix<double> > independent_subgroups_;
};

class vpgl_nitf_RSM_camera_extractor
{
 public:
  // possible outcomes depending on header layout
  enum tre_status{IMAGE_SUBHEADER_TREs_ONLY=0, IMAGE_SUBHEADER_TREs_RSM_TREs,
                  IMAGE_SUBHEADER_TREs_RSM_TREs_OVRFL, INVALID};
   vpgl_nitf_RSM_camera_extractor() = default;
  // path to a NITFV2.1 image file
  vpgl_nitf_RSM_camera_extractor(std::string const& nitf_image_path,
                            bool verbose = false);

  //: Construct from a nitf image pointer
  vpgl_nitf_RSM_camera_extractor(vil_nitf2_image* nift_image,
                            bool verbose = false);

  // image name 
  std::string image_id(size_t image_subheader_index) 
    {
      if(rsm_meta_.count(image_subheader_index)>0){
        rsm_metadata rm = rsm_meta_[image_subheader_index];
        if(rm.image_name_valid)
          return rm.image_name_;
      }
      return "";
    }


  //: number of image subheaders that contain RSM camera TREs
  //  a return of 0 indicates no RSM data
  size_t nitf_header_contains_RSM_tres() const {
      size_t n_RSM = 0;
      for (auto itr = nitf_status_.begin(); itr != nitf_status_.end(); ++itr)
          if (itr->second != INVALID && itr->second != IMAGE_SUBHEADER_TREs_ONLY)
            n_RSM++;
      return n_RSM;
  }

  //: read NITF2.1 tagged record extensions (tres) from header
  // and output a text file of tres present in header
  bool scan_for_RSM_data(bool verbose);
  // text stream of records found
  std::stringstream tre_stream() {
      std::stringstream s(ss_.str());
      return s;
  }
  //: set params for image subheaders that have RSM data
  bool set_RSM_camera_params();

  //: extract adjustable parameter covariance
  bool extract_adjustable_parameter_data();

  bool process(bool verbose){
    if(!scan_for_RSM_data(verbose))
      return false;

    if(!set_RSM_camera_params())
      return false;

    if(!extract_adjustable_parameter_data())
      return false;

    return true;
  }

  //: return the RSM camera associated with the subheader index
  bool RSM_camera(vpgl_RSM_camera<double>& rsm_cam, size_t image_subheader_index = 0 ){
    if(RSM_cams_.count(image_subheader_index) > 0){
      rsm_cam = RSM_cams_[image_subheader_index];
      return true;
    }
    std::cout << "image_subheader index " << image_subheader_index
              << " has no RSM metadata" << std::endl;
    return false;
  }
  // in case of multiple image subheaders the first index associated with
  // a RSM camera definition
  int first_index_with_RSM() {
      if (RSM_cams_.size() == 0)
          return -1;
      auto iter = RSM_cams_.begin();
      return iter->first;
  }
  //: extracted metadata contained in the image subheader including RSM-related info
  // default header index 0
  bool meta( rsm_metadata& rsm_meta,
            ichipb_data& ichipb, size_t image_subheader_index = 0)  {
    if(rsm_meta_.count(image_subheader_index)>0){
      rsm_meta = rsm_meta_[image_subheader_index];
      ichipb = ichipb_data_[image_subheader_index];
      return true;
    }
    std::cout << "image_subheader index " << image_subheader_index <<
      "has no RSM metadata" << std::endl;
    return false;
  }
  bool adjustable_parameter_data(adjustable_parameter_covar& covar_data, size_t image_subheader_index = 0){
    covar_data = adj_param_data_[image_subheader_index];
    if(!covar_data.defined_)
      return false;
    return true;
  }
  // describe the layout of the file header in terms of number of image
  // subheaders and overflow conditions
  void print_file_header_summary();

 private:
  // internal functions
  void ASC_int(std::string str, int& ival){
    std::stringstream ss(str);
    ss >> ival;
  }
  void ASC_double(std::string str, double& dval){
    std::stringstream ss(str);
    ss >> dval;
  }
  // parse the image header tres for required information
  bool determine_header_status(vil_nitf2_image_subheader* header_ptr, size_t header_idx,
      bool& header_has_tres, bool& header_has_RSM, int& ixofl);

  // parse the overflow area for possible RSM information
  bool determine_overflow_status(vil_nitf2_image* nitf_image, size_t header_idx, int ixsofl,
      bool& overflow_has_RSM);

  // determine the layout of information in the file header
  bool init(vil_nitf2_image* nitf_image, bool verbose);

  // extract numerical geographic locations of image corners from
  // the concatenated string representation
  bool process_igeolo(size_t image_subheader_index);

  // data members
  // tres extracted from the image subheader
  std::map<size_t, vil_nitf2_tagged_record_sequence> hdr_ixshd_tres_;

  // tres extracted from the overflow area
  std::map<size_t, vil_nitf2_tagged_record_sequence> ovfl_ixshd_tres_;

  bool RSM_defined_ = false;

  // overall layout status of the file header
  std::map<size_t, tre_status> nitf_status_;

  // useful info from the tres
  std::map<size_t,rsm_metadata> rsm_meta_;

  // cropped image offset
  std::map<size_t,ichipb_data> ichipb_data_;

  // RSM cameras associated with potentially multiple image subheaders
  std::map<size_t, vpgl_RSM_camera<double> > RSM_cams_;

  std::map<size_t, adjustable_parameter_covar> adj_param_data_;

  // Flags indicating if various required TRE groups are present
  bool RSMIDA = false, RSMPIA = false, RSMPCA = false, RSMECA = false, RSMECB = false;

  // presence is checked by attemping to read the EDITION field (40 bytes)
  bool  RSMGIA = false , RSMDCA = false, RSMDCB = false; 
  bool RSMAPA = false, RSMAPB = false, RSMGGA = false;
  int manditory_PCA_row_ = -1; 
  int manditory_PCA_col_ = -1;
  std::stringstream ss_;
};

#endif // vpgl_nitf_RSM_camera_extractor_h_
