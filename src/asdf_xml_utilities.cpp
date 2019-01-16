// TODO - put some utilities in here :-)

#include <iostream>
#include "boost/property_tree/xml_parser.hpp"
#include "boost/foreach.hpp"
#include "asdf_xml_utilities.h"

namespace pt = boost::property_tree;

int parse_quakeml_data(const char* quakeml_string, struct quakeml& qml_data)
{
  qml_data.nevents = -1;
  std::istringstream qml_is(quakeml_string);
  pt::ptree qml_tree;
  try {
    pt::read_xml(qml_is, qml_tree);
    // printf("QuakeML latitude:\n%s\n", qml_tree.get<std::string>("q:quakeml.eventParameters.event.origin.time.value").c_str());
    // double longitude = qml_tree.get<double>("q:quakeml.eventParameters.event.origin.latitude.value");
    // printf("QuakeML latitude:\n%1.16f\n", longitude );
  } 
  catch (...)
  {
    std::cerr << "Warning, QuakeML parsing from ASDF file failed!" 
      << std::endl;
    return -1;
  }

  // Get the eventParameters publicID
  std::string prop = "q:quakeml.eventParameters.<xmlattr>.publicID";
  try {
    std::string eP_id = qml_tree.get<std::string>(prop);
    qml_data.eventParameters_publicID = eP_id;
    // std::cout << "QuakeML eP_ID: " << eP_id << std::endl;
  } 
  catch (...)
  {
    std::cerr << "Warning, QuakeML section of ASDF file missing: "
      << prop << std::endl;
    return -1;
  }

  // Get the events and their publicID attributes
  try
  {
    prop = "q:quakeml.eventParameters.event.<xmlattr>";
    int nevents = 0;
    BOOST_FOREACH(pt::ptree::value_type &v, qml_tree.get_child(prop))
    {
      std::string attr = v.first.data();
      if (attr.compare("publicID"))
        continue;
      // Parse the eventID
      std::string ev_id = v.second.data();
      qml_data.event_publicID.push_back(ev_id);
      // std::cout << "  QuakeML event id: " << ev_id << std::endl;
      nevents++;
    }
    std::cout << "Read QuakeML: " << nevents << " events" << std::endl;
    qml_data.nevents = nevents;
  } 
  catch (...)
  {
    std::cerr << "Warning, could not read QuakeML events section of ASDF"
      << std::endl;
    return -1;
  }

  return qml_data.nevents;
}


int parse_stationxml_data(const char* staxml_char, 
    struct stationxml& staxml_data)
{
  std::istringstream sml_is(staxml_char);
  pt::ptree sml_tree;
  try {
    pt::read_xml(sml_is, sml_tree);
    // printf("QuakeML latitude:\n%s\n", qml_tree.get<std::string>("q:quakeml.eventParameters.event.origin.time.value").c_str());
    std::string prop("FDSNStationXML.Network.<xmlattr>.code");
    staxml_data.name = sml_tree.get<std::string>(prop);
    staxml_data.name.append(".");
    prop = "FDSNStationXML.Network.Station.<xmlattr>.code";
    staxml_data.name.append(sml_tree.get<std::string>(prop));
    prop = "FDSNStationXML.Network.Station.Latitude";
    staxml_data.latitude = sml_tree.get<double>(prop);
    prop ="FDSNStationXML.Network.Station.Longitude";
    staxml_data.longitude = sml_tree.get<double>(prop);
    prop ="FDSNStationXML.Network.Station.Elevation";
    staxml_data.elevation = sml_tree.get<double>(prop);
  } 
  catch (...)
  {
    std::cerr << "Warning, StationXML parsing from ASDF file failed!" 
      << std::endl;
    return -1;
  }

  return 0;
}
