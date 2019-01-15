/**
 * @file read_ASDF.cpp
 */

#include <hdf5.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "boost/property_tree/xml_parser.hpp"
#include "ASDF_common.h"
#include "ASDF_init.h"
#include "ASDF_read.h"
#include "ASDF_write.h"

#define FAKE_NUM_SAMPLES 30
#define MAX_STRING_LENGTH 256

namespace pt = boost::property_tree;

int main(int argc, char *argv[]) {
  // int i;

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /******************************************************
   *   Begin Fake Data                                  *
   ******************************************************/
  // std::string filename = "new_file.h5";
  std::string filename = "preprocessed.h5";
  /*
  char *event_name = "event0123456789";
  char *station_name = "AF.CVNA";
  char *station_xml = "station_xml_string";
  char *provenance_string = "provenance_string";
  char *quakeml_string = "quakemlstringstring";
  */
    /*
      "quakemlstring = '<quakeml>\\n<event unique_id=\"EV_01\">\\n<location"
      "main=\"true\" unique_id=\"LOC_01\" analysis-type=\"M\">\\n<origin-date"
      "timezone=\"00:00\">\\n<year\\>2004</year>\\n<month\\>09</month\\n"
      "<day\\>28</day>\\n<hour\\>17</hour>\\n<minute\\>15</minute>\\n"
      "<seconds\\>24.0</seconds>\\n</origin-date>\\n<latitude error=\"0\">"
      "35.8</latitude>\\n<longitude error=\"0\">-120.4</longitude>\\n<depth"
      "unit=\"km\" error=\"0\">7</depth>\\n<magnitude unit=\"M\""
      "error=\"0\">6.0</magnitude>\\n<region>CENTRAL"
      "CALIFORNIA</region>\\n<author>SPECFEM3D_GLOBE</author>\\n</location>\\n"
      "</event>\\n</quakeml>";
    */

/*
void ASDFIO::readASDF( vector<TimeSeries*>& a_TimeSeries)
{
#ifdef USE_ASDF
*/
  printf("Reading ASDF file: %s\n", filename.c_str());

  ASDF_initialize_hdf5();

/*
  stringstream filePrefix;
  string path = m_ew->getPath();
//building the file name...
  if( path != "." )
    filePrefix << path;
  filePrefix << m_fileName;
*/

  hid_t file_id;
  file_id = ASDF_open_read_only(filename.c_str(), MPI_COMM_WORLD);
  if (file_id < 0)
  {
    printf("Could not open ASDF hdf5 file: %s\n", filename.c_str());
    MPI_Abort(MPI_COMM_WORLD, file_id);
  }

  char* quakeml = NULL;
  herr_t err = ASDF_read_quakeml(file_id, &quakeml);
  if (err < 0)
  {
    printf("Could not read QuakeML in ASDF hdf5 file: %s\n", filename.c_str());
    MPI_Abort(MPI_COMM_WORLD, file_id);
  }
  printf("QuakeML contents:\n%s\n", quakeml);
  pt::ptree qml_tree;
  std::string qml(quakeml);
  std::istringstream qml_is(qml);
  pt::read_xml(qml_is, qml_tree);
  // printf("QuakeML latitude:\n%s\n", qml_tree.get<std::string>("q:quakeml.eventParameters.event.origin.time.value").c_str());
  double longitude = qml_tree.get<double>("q:quakeml.eventParameters.event.origin.latitude.value");
  printf("QuakeML latitude:\n%1.16f\n", longitude );
  // printf("QuakeML latitude:\n%e\n", qml_tree.get("q:quakeml.eventParameters.event.origin.latitude",0.0));

  hid_t waveforms_grp = ASDF_open_waveforms_group(file_id);
  hsize_t num_stn;
  err = H5Gget_num_objs(waveforms_grp, &num_stn);
  printf("\nWaveform stations: %d\n", (int) num_stn);
  for (int stn=0; stn < num_stn; stn++)
  {
    char buff[1000];
    err = H5Gget_objname_by_idx(waveforms_grp, (hsize_t) stn, buff, (size_t) 1000);
    printf("  station %d: %s\n", stn, buff);

    char* stationxml = NULL;
    err = ASDF_read_stationxml(waveforms_grp, std::string(buff).c_str(), 
        &stationxml);

    printf("    stationxml:\n%s\n", stationxml);
    pt::ptree sml_tree;
    std::string sml(stationxml);
    std::istringstream sml_is(sml);
    pt::read_xml(sml_is, sml_tree);
    double longitude = sml_tree.get<double>("FDSNStationXML.Network.Station.Latitude");
    printf("    station latitude:\n%1.16f\n", longitude );
    // printf("    station latitude:\n%s\n", qml_tree.get<std::string>("FDSNStationXML.Source").c_str());
    free(stationxml);
  }
/*

  // Only support 3 xyz or vel components output for now
  const int num_waveforms=3;
  // Allocate char* buffers
  const int buff_size=100;
  char **waveform_names = (char **) malloc(num_waveforms*sizeof(char *));
  char *event_name = (char *) malloc(buff_size*sizeof(char));
  sprintf(event_name, "UNKNOWN");
  for (int i=0; i < num_waveforms; ++i)
    waveform_names[i] = (char *) malloc(buff_size*sizeof(char));

  // Allocate a temp that can hold our biggest time series data
  int max_size = 0;
  for (int ts=0; ts < a_TimeSeries.size(); ts++)
    max_size = max(max_size, a_TimeSeries[ts]->getNsteps());
  float *data = (float *) malloc(max_size*sizeof(float));

  // For each time series, create the station group
  for (int ts=0; ts < a_TimeSeries.size(); ts++)
  {
    // Skip if it's not our point
    bool myPoint = a_TimeSeries[ts]->myPoint();
    if (!myPoint)
      continue;

    string station_name = a_TimeSeries[ts]->getStationName();
    TimeSeries::receiverMode mode = a_TimeSeries[ts]->getMode();
    if (!(mode == TimeSeries::Displacement) && !(mode == TimeSeries::Velocity))
    {
      cout << "ASDF hdf5 outputs only displacement or velocity ... skipping: " 
        << station_name << endl;
      continue;
    }

    hid_t station_grp = 
      ASDF_open_stations_group(waveforms_grp, station_name.c_str());
    if (station_grp < 0)
    {
      cout << "Could not find ASDF hdf5 station group: " << station_name
        << " ... skipping." << endl;
      continue;
    }

    hid_t data_id[num_waveforms];
    int nsamples=a_TimeSeries[ts]->getNsteps();
    int start_time=a_TimeSeries[ts]->get_shift(); // TODO - right value?
    double sampling_rate=m_ew->getTimeStep(); // TODO - right value?
    char disp_names[3][2] = {"x","y","z"}; // 2 for ?\0
    char vel_names[3][3] = {"ux","uy","uz"}; // 3 for ??\0
    for (int i=0; i < num_waveforms; ++i)
    {
      if (mode == TimeSeries::Displacement)
        sprintf(waveform_names[i], "displacement_%s", disp_names[i]);
      else
        sprintf(waveform_names[i], "velocity_%s", vel_names[i]);
      data_id[i] = ASDF_open_waveform(station_grp, waveform_names[i]);
    }

    assert(nsamples);
    int chunk_size = 10; // TODO - what's best read chunk?
    int num_chunk = (nsamples - 1) / chunk_size + 1;

    float_sw4** ts_data = a_TimeSeries[ts]->getRecordingArray();
    for (int i = 0; i < num_waveforms; ++i)
    {
      // ASDF_write_partial_waveform(data_id[i], data, offset, samples);
      H5Dread(data_id[i], H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      for (int j=0; j < nsamples; j++)
          ts_data[i][j] = data[j];
    }
#if 0
#endif

    for (int i = 0; i < num_waveforms; ++i)
      H5Dclose(data_id[i]);

    ASDF_close_group(station_grp);
    // For testing / comparing ASDF and USGS output
    a_TimeSeries[ts]->writeFile();
  }

  // Free mallocs
  for (int i = 0; i < num_waveforms; ++i)
    free(waveform_names[i]);
  free(data);
  free(waveform_names);
  free(event_name);
*/
  free(quakeml);

  // Close everything
  ASDF_close_group(waveforms_grp);
  H5Fclose(file_id);

  ASDF_finalize_hdf5();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
