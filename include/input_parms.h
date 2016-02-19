#ifndef INPUT_PARMS_H_
#define INPUT_PARMS_H_

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include "macros.h"

namespace LapH {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class input_parameter{
public:
  // ---------------------- METHODS ------------------------------------------
  // standard constructor, does not do anything
  input_parameter() {}; 
  // standard destructor, does not do anything
  ~input_parameter() {}; 

public:
  // --------------------- INTERFACE -----------------------------------------
  void parse_input_file(int argc, char *argv[]);
  void print_options();
  // -------------------- END INTERFACE --------------------------------------

private:
  //  DISALLOW_COPY_AND_ASSIGN(input_parameter); // no copy nor assignment
  // checks the input parameters for consistency
  void check_input_parameters();
  // creates the filenames for the perambulator and the random vectors and
  // checks if these files already exist. If they exist the program is
  // terminated to avoid overwriting these files!
  void check_and_create_filenames();

public:
  // -------------------------- DATA -----------------------------------------
  size_t dilution_size_so[3];      // SOURCE dilution sizes and types
  std::string dilution_type_so[3]; // 0 = time, 1 = eigevector space, 2 = Dirac
  size_t dilution_size_si[4];      // SINK dilution sizes and types
  std::string dilution_type_si[4]; // 0 = time, 1 = space, 2 = Dirac, 3 = color
                                   // NOTE: Please use capital Latters I, B, F!!
                                   // All of these numbers are NOT the blocksize
                                   // or interlace sizes, but the number of 
                                   // inversions in the corresponding space 
  size_t config;                   // configuration number
  size_t nb_config;                // total number of configs to be processed
  size_t delta_config;             // distance between configurations
  size_t Ls;                       // spatial extend
  size_t Lt;                       // temporal extend
  size_t nb_ev;                    // number of eigenvectors
  size_t nb_rnd;                   // random vector id of vector in process
  size_t verbose;                  // displaying more informations

  size_t peram_gen_omp_num_threads;// number of OpenMP threads

  int seed;                        // seed for random vector generation

  std::string quarktype;           // quark type: u,d,s,c -> for naming outfiles
  std::string outpath;             // path to write everything
  std::string inpath_ev;           // path to read eigenvectors
  std::string peram_file_name;     // perambulator file name
  std::string rnd_vec_file_name;   // random vektor file name
  // ------------------------- END DATA ---------------------------------------

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ends here

#endif /* INPUT_PARMS_H_ */

