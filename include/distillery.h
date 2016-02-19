#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <complex>
#include <iostream>
#include <fstream>
#include <string>

#include "mpi.h"

#include <Eigen/Dense>
#include <Eigen/Core>

#include "macros.h"
#include "ranlxs.h"
#include "input_parms.h"

namespace LapH {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class distillery{
  // ------------------------ METHODS ------------------------------------------
public:
  // standard constructor - doesn't do much -> Work is done in this->initialise
  distillery(){
    V = NULL;
    perambulator  = NULL;
    random_vector = NULL;
    memory_allocated = false;
  };
  // standard deconstructor
  ~distillery(){ 
    if(memory_allocated == true)
      this->clean();
  };

  // ---------------------------------------------------------------------------
  // ------------------------- INTERFACE ---------------------------------------
  // Initialisation:
  // allocation of memory, reading eigenvectors, generating randomvector
  // input: in_param -> all necessary input parameters (see struct definition)
  void initialise(const input_parameter& in_param);
  // cleaning everything up and freeing all memory
  void clean();
  // setting everything to zero again, and reding new eigenvectors
  // only difference to initialise is that no new memory is allocated
  // input: in_param -> all necessary input parameters (see struct definition)
  void reset_all(const input_parameter& in_param);
  // setting the perambulators to zero again, e.g. for a new rnd vector run
  // and generating a new random vector
  // input: seed -> new seed to generate new random vector
  //        random_vector_number -> just the id of random vector in progress
  void reset_perambulator_and_randomvector(const int seed, 
                                           const size_t random_vector_number);
  // creates the source for inversions -- param.dilution_type_so[2] sources are
  // created simultaneously. This happens because the sources in Dirac space
  // are very similar and it is faster to create these sources at once
  // input:   dil_t -> counts the number of inversions in time
  //          dil_e -> counts the number of inversions in eigenvector space
  // output:  source -> complex vector to store the sources 
  void create_source(const size_t dil_t, const size_t dil_e, 
                     std::complex<double>* source);
  // multiplies a propagator with V^dagger and stores it in one column of
  // of the perambulator. The column number is computed with dil_t, dil_e,
  // and dil_d
  // input:   dil_t -> counts the number of inversions in time
  //          dil_e -> counts the number of inversions in eigenvector space
  //          dil_d -> counts the number of inversions in Dirac space
  //          propagator -> the propagator contains the inverted Dirac
  //                        operator
  void add_to_perambulator(const size_t dil_t, 
                           const size_t dil_e, 
                           const size_t dil_d, 
                           const std::complex<double>* const propagator);
  // writes the perambulator to disk
  // the indices from slow to fast: t -> Dirac -> eigenvectors :: This ordering
  // is important because it speeds up the multiplication with V^dagger which
  // is blockdiagonal in Dirac space
  void write_perambulator_to_disk();
  // ----------------------- END INTERFACE -------------------------------------
  // ---------------------------------------------------------------------------

private:
  DISALLOW_COPY_AND_ASSIGN(distillery); // no copy nor assignment
  void read_eigenvectors();
  void set_random_vector();
  void read_random_vector();
  void write_random_vector_to_disk();
  // -------------------------- DATA -------------------------------------------
  Eigen::MatrixXcd* V;             // memory eigensystem
  Eigen::MatrixXcd* perambulator;  // memory perambulator
  Eigen::VectorXcd* random_vector; // memory random vector
  input_parameter param;           // all necessary input parameters

  // some flags
  bool memory_allocated;           // is all memory allocated?

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}

#endif /* BASICOPERATOR_H_ */
