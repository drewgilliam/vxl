//-*- c++ -*-------------------------------------------------------------------
#ifdef __GNUC__
#pragma implementation
#endif
#include "vbl_bool_ostream.h"


ostream& operator<<(ostream& s, const vbl_bool_ostream::on_off& proxy) {
  if (*(proxy.truth)) 
    s << "on";
  else 
    s << "off";
  return s;
}
 
ostream& operator<<(ostream& s, const vbl_bool_ostream::high_low& proxy) {
  if (*(proxy.truth)) 
    s << "high";
  else 
    s << "low";
  return s;
}
 
ostream& operator<<(ostream& s, const vbl_bool_ostream::true_false& proxy) {
  if (*(proxy.truth)) 
    s << "true";
  else 
    s << "false";
  return s;
}
 
