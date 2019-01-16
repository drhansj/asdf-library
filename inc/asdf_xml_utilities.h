/**
 * @file asdf_xml_utilities.h
 * @brief
 * @author Hans Johansen
 */

#ifndef _ASDF_XML_UTILITIES_H_
#define _ASDF_XML_UTILITIES_H_

/**
 * @brief A structure to hold minimal QuakeML data
 */
struct quakeml {
  int nevents;
  std::string eventParameters_publicID;
  std::vector<std::string> event_publicID;
};  // struct quakeml

int parse_quakeml_data(const char* quakeml_string, 
    struct quakeml& quakeml_data);

/**
 * @brief A structure to hold minimal StationXML data
 */
struct stationxml {
  std::string name; // NET.STA name
  float latitude;
  float longitude;
  float elevation;
};  // struct stationxml

int parse_stationxml_data(const char* staxml_string, 
    struct stationxml& staxml_data);

#endif  // _ASDF_XML_UTILITIES_H_
