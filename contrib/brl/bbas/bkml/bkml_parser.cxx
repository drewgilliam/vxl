//:
// \file
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstddef>
#include "bkml_parser.h"
//
// \brief Parses the kml configuration file for bwm tool.
// \verbatim
//  Modifications
//   2012-09-10 Yi Dong - Modified to parser the polygon and path(LineString) coordinates stored in kml
// \endverbatim
//
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif

// --------------
// --- PARSER ---
// --------------
template <typename T>
void convert(const char* t, T& d)
{
  std::stringstream strm(t);
  strm >> d;
}

bkml_parser::bkml_parser()
{
  init_params();
}

void bkml_parser::init_params()
{
  heading_ = 0.0;
  tilt_ = 0.0;
  roll_ = 0.0;
  right_fov_ = 0.0;
  top_fov_ = 0.0;
  near_ = 0.0;
  heading_dev_ = 0.0;
  tilt_dev_ = 0.0;
  roll_dev_ = 0.0;
  right_fov_dev_ = 0.0;
  top_fov_dev_ = 0.0;
  tags_.clear();
  data_.clear();
}

// parse from filename
bool
bkml_parser::parseFilename(const std::string& filename)
{
  // open file
  std::FILE* file_ptr = std::fopen(filename.c_str(), "r");
  if (!file_ptr) {
    std::cerr << "ERROR opening KML " << filename << std::endl;
    return false;
  }

  // parse file via parent function
  if (!this->parseFile(file_ptr)) {
    std::cerr << "ERROR for KML " << filename << std::endl
              << XML_ErrorString(this->XML_GetErrorCode()) << " at line "
              << this->XML_GetCurrentLineNumber() << std::endl;
    return false;
  }

  // cleanup
  return true;
}


void
bkml_parser ::cdataHandler(const std::string&  /*name*/, const std::string&  /*data*/)
{
}

void
bkml_parser::handleAtts(const XML_Char** /*atts*/)
{
}

// begin element parsing
void
bkml_parser::startElement(const char* name, const char**  /*atts*/)
{
  // record element list of tags
  this->tags_.push_back(name);
}

// character data for a single tag may be split among multiple charData calls
// https://libexpat.github.io/doc/common-pitfalls/#split-character-data
// This function acculumates character data, later parsed in endElement
void bkml_parser::charData(const XML_Char* s, int len)
{
  const int leadingSpace = skipWhiteSpace(s);
  if (len==0 || len<=leadingSpace)
    return;  // called with whitespace between elements

  // accumulate character data
  for (int i = 0; i<len; ++i)
    this->data_ += s[i];
}

// end element parsing (wrapper for processElement)
void
bkml_parser::endElement(const char* name)
{
  // name as std::string
  std::string tag(name);

  // verbose output
  // std::cout << "  tag list: ";
  // for (auto t : this->tags_)
  //   std::cout << t << ",";
  // std::cout << std::endl;
  //
  // std::cout << "CURRENT TAG: " << tag << std::endl;
  // std::cout << "DATA:" << this->data_ << std::endl;

  // check for malformed KML
  if (this->tags_.back() != tag) {
    std::cerr << "ERROR: KML may be malformed, unexpected tag end </" << tag << ">" << std::endl;
    this->data_.clear();
    return;
  }

  // process non-empty data
  else if (!this->data_.empty()) {
    this->processElement();
  }

  // clear tag & character data
  this->tags_.pop_back();
  this->data_.clear();
}

// process tag/data pair
void
bkml_parser::processElement()
{
  // current tag
  std::string tag = this->tags_.back();

  // istringstream for parsing
  std::istringstream iss(this->data_);

  // simple parsing
  if (tag == KML_LON_TAG) {
    iss >> longitude_;
  } else if (tag == KML_LAT_TAG) {
    iss >> latitude_;
  } else if (tag == KML_PLACEMARK_NAME_TAG) {
    current_name_ = iss.str();
  } else if (tag == KML_ALT_TAG) {
    iss >> altitude_;
  } else if (tag == KML_HEAD_TAG) {
    iss >> heading_;
  } else if (tag == KML_HEAD_DEV_TAG) {
    iss >> heading_dev_;
  } else if (tag == KML_TILT_TAG) {
    iss >> tilt_;
  } else if (tag == KML_TILT_DEV_TAG) {
    iss >> tilt_dev_;
  } else if (tag == KML_ROLL_TAG) {
    iss >> roll_;
  } else if (tag == KML_ROLL_DEV_TAG) {
    iss >> roll_dev_;
  } else if (tag == KML_RFOV_TAG) {
    iss >> right_fov_;
  } else if (tag == KML_RFOV_DEV_TAG) {
    iss >> right_fov_dev_;
  } else if (tag == KML_TFOV_TAG) {
    iss >> top_fov_;
  } else if (tag == KML_TFOV_DEV_TAG) {
    iss >> top_fov_dev_;
  } else if (tag == KML_NEAR_TAG) {
    iss >> near_;
  }

  // coordinate parsing
  else if (tag == KML_CORD_TAG) {

    // tag & vertices
    std::string coord_tag;
    std::vector<vgl_point_3d<double> > poly;

    // discover coordinate tag
    // (iterate backwards through tags_  to find first recognized identifier)
    for (auto it = this->tags_.rbegin(); it != this->tags_.rend(); ++it ) {

      // check current tag against known tags
      if (*it == KML_POLYOB_TAG)
        coord_tag = KML_POLYOB_TAG;
      else if (*it == KML_POLYIB_TAG)
        coord_tag = KML_POLYIB_TAG;
      else if (*it == KML_LINE_TAG)
        coord_tag = KML_LINE_TAG;
      else if (*it == KML_POINT_TAG )
        coord_tag = KML_POINT_TAG;

      // if known tag, we're done!
      if (!coord_tag.empty())
        break;
    }

    // quit on invalid coord_tag
    if (coord_tag.empty()) {
      std::cerr << "Unrecognized coordinate tag" << std::endl;
      return;
    }

    // parse coordinates
    // --coordinate elements separated by commas
    // --coordinate tuple/triple separated by whitespace
    // --for example "52.3,83.2,0.0  52.4,83.4,0.0"
    std::vector<double> coord;
    double val;

    while (iss >> val) {

      // comma separated coordiante elements
      coord.push_back(val);
      if (iss.peek() == ',') {
        iss.ignore(1);
        continue;
      }

      // save coordinate to poly
      coord.resize(3,0.0);
      poly.push_back(vgl_point_3d<double>(coord[0],coord[1],coord[2]));

      // clear for next coordinate
      coord.clear();
    }

    // check for valid coordinates
    if (poly.empty()) {
      std::cerr << "Empty coordinates for <" << coord_tag << ">" << std::endl;
      return;
    }

    // VXL polygons are implictly closed
    // thus remove the last element if equal to first element
    if (poly.front() == poly.back()) {
      poly.pop_back();
    }

    // save coordinates
    if (coord_tag == KML_POLYOB_TAG)
      polyouter_.push_back(poly);
    else if (coord_tag == KML_POLYIB_TAG)
      polyinner_.push_back(poly);
    else if (coord_tag == KML_LINE_TAG)
      linecord_.push_back(poly);
    else if (coord_tag == KML_POINT_TAG )
      points_.push_back(poly[0]);

  } // end coordinate parsing

  // unrecognized tag
  else {
    // std::cerr << "Unrecognized tag <" << tag << ">" << std::endl;
  }

}

void bkml_parser::trim_string(std::string& s)
{
  int i = (int)s.find_first_not_of(' ');
  int j = (int)s.find_last_not_of(' ');
  std::string t = s.substr(i,j-i+1);
  s = t;
}

// parse points given a kml filename
std::vector<vgl_point_3d<double> >
bkml_parser::parse_points(const std::string& filename)
{
  // initialize outputs
  std::vector<vgl_point_3d<double> > out;
  out.clear();

  // parse file
  bkml_parser parser;
  if (!parser.parseFilename(filename))
    return out;

  // cleanup
  return parser.points_;
}

// parse polygon with outer boundary given a kml filename
vgl_polygon<double>
bkml_parser::parse_polygon(const std::string& filename)
{
  // initialize outputs
  vgl_polygon<double> out;
  out.clear();

  // parse file
  bkml_parser parser;
  if (!parser.parseFilename(filename))
    return out;

  // check for empty polygon
  if (parser.polyouter_.empty()) {
    std::cerr << "input kml has no polygon outerboundary, return an empty polygon" << '\n';
    return out;
  }

  // load the outer boundary
  for (auto & sh_idx : parser.polyouter_) {
    out.new_sheet();
    auto n_points = (unsigned)sh_idx.size();
    for (unsigned pt_idx = 0; pt_idx < n_points; pt_idx++) {
      out.push_back(sh_idx[pt_idx].x(), sh_idx[pt_idx].y());
    }
  }

  // cleanup
  return out;
}

// parse polygon with outer & inner boundary given a kml filename
vgl_polygon<double>
bkml_parser::parse_polygon_with_inner(const std::string& filename, vgl_polygon<double>& outer, vgl_polygon<double>& inner,
                                      unsigned& n_out, unsigned& n_in)
{
  // initialize outputs
  vgl_polygon<double> out;
  out.clear();
  outer.clear();
  inner.clear();

  // parse file
  bkml_parser parser;
  if (!parser.parseFilename(filename))
    return out;

  // check for empty polygon
  if (parser.polyouter_.empty()) {
    std::cerr << "input kml has no polygon outerboundary, return an empty polygon" << '\n';
    return out;
  }

  // create polygon from parser
  n_out = (unsigned)parser.polyouter_.size();
  n_in = (unsigned)parser.polyinner_.size();

  // load the outer boundary
  for (unsigned sh_idx = 0; sh_idx < n_out; sh_idx++) {
    out.new_sheet();
    outer.new_sheet();
    auto n_points = (unsigned)parser.polyouter_[sh_idx].size();
    for (unsigned pt_idx = 0; pt_idx < n_points; pt_idx++) {
      out.push_back(parser.polyouter_[sh_idx][pt_idx].x(), parser.polyouter_[sh_idx][pt_idx].y());
      outer.push_back(parser.polyouter_[sh_idx][pt_idx].x(), parser.polyouter_[sh_idx][pt_idx].y());
    }
  }

  // load the inner boundary
  for (unsigned sh_idx = 0; sh_idx < n_in; sh_idx++) {
    out.new_sheet();
    inner.new_sheet();
    auto n_points = (unsigned)parser.polyinner_[sh_idx].size();
    for (unsigned pt_idx = 0; pt_idx < n_points; pt_idx++) {
      out.push_back(parser.polyinner_[sh_idx][pt_idx].x(), parser.polyinner_[sh_idx][pt_idx].y());
      inner.push_back(parser.polyinner_[sh_idx][pt_idx].x(), parser.polyinner_[sh_idx][pt_idx].y());
    }
  }

  // cleanup
  return out;
}

// parse single location given a kml filename
bool
bkml_parser::parse_location_from_kml(const std::string& filename, double& lat, double& lon)
{
  // initialize outputs
  lat = 0.0;
  lon = 0.0;

  // parse file
  bkml_parser parser;
  if (!parser.parseFilename(filename))
    return false;

  // cleanup
  lat = parser.latitude_;
  lon = parser.longitude_;
  return true;
}
