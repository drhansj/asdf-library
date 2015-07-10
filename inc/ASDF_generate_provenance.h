/******************************************************************************
 * Copyright 2015 ASDF developers
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/
/** 
 * @file ASDF_generate_provenance.h
 * @brief 
 * @author Matthieu Lefebvre
 */

#ifndef _ASDF_GENERATE_PROVENANCE_H_
#define _ASDF_GENERATE_PROVENANCE_H_

/**
 *  @brief Generate a provenance entity with prov_id and prov_label
 *         for the software used in the waveform simulation.
 *
 *  @param software_name Name of the solver, for instance "SPECFEM3D_GLOBE"
 *  @param software_version Version of the solver, for instance "7.0"
 *  @param software_website Website that hosts the solver
 *  @param prov_label Description of the provenance.
 *  @param prov_id Unique id for this provenance item.
 *  @param prov Output, pointer to the generated C string.
 */
void ASDF_generate_software_agent(const char *software_name,
                                  const char *software_version,
                                  const char *software_website,
                                  const char *prov_label,
                                  const char *prov_id,
                                  char **prov);

/**
 *  @brief Generate a provenance entity with prov_id and prov_label
 *         for the waveform trace generated from the solver
 *
 *  @param prov_label Description of the provenance.
 *  @param prov_id Unique id for this provenance item.
 *  @param prov Output, pointer to the generated C string.
 */
void ASDF_generate_trace_entity(const char *prov_label,
                                const char *prov_id,
                                char **prov);

/**
 *  @brief Generate a provenance activity with prov_id and prov_label
 *         for the simulation
 *
 *  @param startTime start time of the simulation, for instance 2014-02-02T12:15:03
 *  @param endTime end time of the simulation, for instance 2014-02-02T14:07:13
 *  @param prov_label Description of the provenance.
 *  @param prov_id Unique id for this provenance item.
 *  @param prov Output, pointer to the generated C string.
 */
void ASDF_generate_simulation_activity(const char *startTime,
                                       const char *endTime,
                                       const char *prov_label,
                                       const char *prov_id,
                                       char **prov);


/**
 *  @brief Generate a provenance association using the unique prov_ids such that 
 *         prov_id1 is associated with prov_id2
 *
 *  @param prov_id1 Unique id for the associated provenance item.
 *  @param prov_id2 Unique id for the associated provenance item.
 */
void ASDF_generate_association(const char *prov_id1,
                               const char *prov_id2,
                               char **prov);

/**
 *  @brief Generate a provenance usage using the unique prov_ids such that 
 *         prov_id1 used prov_id2
 *
 *  @param prov_id1 Unique id for the used provenance item.
 *  @param prov_id2 Unique id for the used provenance item.
 */
void ASDF_generate_usage(const char *prov_id1,
                         const char *prov_id2,
                         char **prov);

/**
 *  @brief Generate a provenance generation using the unique prov_ids such that 
 *         prov_id1 was generated by prov_id2
 *
 *  @param prov_id1 Unique id for the generated provenance item.
 *  @param prov_id2 Unique id for the generated provenance item.
 */

void ASDF_generate_wasGeneratedby(const char *prov_id1,
                                  const char *prov_id2,
                                  char **prov);

/** 
 * @brief Generate a provenance entity with id prov_id and label prov_label
 *        from a specfem style parameter file ("Par_file") into a C string prov.
 * 
 * @param filename Name of the "Par_file".
 * @param prov_label Description, for instance "SPECFEM Input Parameters".
 * @param prov_id Unique id for this provenance item.
 * @param prov Output, pointer to the generated C string.
 */
void ASDF_generate_par_file_provenance(const char *filename,
                                       const char *prov_label,
                                       const char *prov_id,
                                       char **prov);

/** 
 * @brief Deallocate the sub_provenance string build by
 *        ASDF_generate_par_file_provenance and other provenance subroutines
 * 
 * @param sub_provenance Pointer to the C string to deallocate.
 */
void ASDF_clean_par_file_provenance(char *sub_provenance);

void ASDF_generate_sf_provenance(const char *filename, const char *startTime, const char *endTime, char **prov);

#endif  /* _ASDF_GENERATE_PROVENANCE_H_ */
