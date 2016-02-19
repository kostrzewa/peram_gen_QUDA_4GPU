#include "input_parms.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// checking that the dilution parameters are correctly set in infile
// input:  type       -> gives the type of the dilution scheme
//         max_size   -> maximal size of space
// output: nb_dil_vec -> sets (in case of full and no dilution) and checks the 
//                       number of dilution vectors
static void check_dilution_input(const std::string type, const size_t max_size, 
                                 size_t& nb_dil_vec){
  // check for type
  if( (type.compare("F") != 0) && 
      (type.compare("I") != 0) && 
      (type.compare("B") != 0) && 
      (type.compare("N") != 0) ) {
        std::cerr << "Inversion type has to be one of \"F\", \"I\"," \
                     " \"B\" or \"N\"." << std::endl;
        std::cerr << "Aborting..." << std::endl;
        exit(0);
  }
  // check and set number of inversions in corresponding space
  // TODO: must be change because it gives a strange result for Sinksmearing!!
  if (type.compare("F") == 0)
     nb_dil_vec = max_size;
  else if (type.compare("N") == 0 )
     nb_dil_vec = 1;
  else {
    if(max_size % nb_dil_vec != 0) {
      std::cerr << "check number of inversions" << std::endl;
      exit(0);
    }
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::check_input_parameters(){

  if(peram_gen_omp_num_threads > 2 || peram_gen_omp_num_threads < 1){
    std::cout << "The number of OMP threads for this version of peram_gen can only be 1 or 2! The selected value was " <<
                 peram_gen_omp_num_threads << std::endl;
    exit(1);
  }

  if(config > 100000) {
    std::cout << "Please check whether configuration number is correct: " 
              << config << std::endl;
    exit(0);
  }
  if(Ls > 500) {
    std::cout << "Please check whether Ls is correct: " <<  Ls << std::endl;
    exit(0);
  }
  if(Lt > 500) {
    std::cout << "Please check whether Lt is correct: " <<  Lt << std::endl;
    exit(0);
  }
  if(nb_ev > 5000) {
    std::cout << "Please check whether number of eigenvectors is correct: " 
              <<  nb_ev << std::endl;
    exit(0);
  }
  if(nb_rnd > 1000) {
    std::cout << "Please check whether starting number of randomvector" \
                 " is correct: " <<  nb_rnd << std::endl;
    exit(0);
  }
  if((quarktype.compare("u") != 0) && 
     (quarktype.compare("d") != 0) && 
     (quarktype.compare("s") != 0) && 
     (quarktype.compare("c") != 0) ) {
       std::cerr << "Quarktype has to be one of \"u\", \"d\", \"s\" or \"c\"." 
                 << std::endl;
       std::cerr << "Aborting..." << std::endl;
       exit(0);
  }
  check_dilution_input(dilution_type_so[0], Lt, dilution_size_so[0]); 
  check_dilution_input(dilution_type_so[1], nb_ev, dilution_size_so[1]); 
  check_dilution_input(dilution_type_so[2], 4, dilution_size_so[2]);
  check_dilution_input(dilution_type_si[0], Lt, dilution_size_si[0]); 
  check_dilution_input(dilution_type_si[1], Ls*Ls*Ls, dilution_size_si[1]); 
  check_dilution_input(dilution_type_si[2], 4, dilution_size_si[2]); 
  check_dilution_input(dilution_type_si[3], 3, dilution_size_si[3]); 
  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::check_and_create_filenames(){

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::parse_input_file(int argc, char *argv[]) {

  int opt = -1;
  int reader = 0;
  char infilename[200];
  char readin[256];
  FILE* infile = NULL;

  // search for command line option and put filename in "infilename"
  for(int i = 0; i < argc; ++i) {
    if(std::strcmp(argv[i], "-LapHsin") == 0) {
      opt = i+1;
      break;
    }
  }
  if(opt < 0) {
    std::cout << "No input file specified, trying infile.in" << std::endl;
    sprintf(infilename, "infile.in");
  } else {
    sprintf(infilename, "%s", argv[opt]);
    std::cout << "Trying input file " << infilename << std::endl;
  }
  // open file for reading
  if ((infile = fopen(infilename, "r")) == NULL ) {
    std::cerr << "Could not open file " << infilename << std::endl;
    std::cerr << "Aborting..." << std::endl;
    exit(-10);
  }

  // scan infile and check arguments ----------------------------------------
  // configs
  peram_gen_omp_num_threads = 2;
  reader += fscanf(infile, "omp_num_threads = %zu\n", &peram_gen_omp_num_threads);
  reader += fscanf(infile, "config = %zu \n", &config);
  reader += fscanf(infile, "total number of configs = %zu \n", &nb_config);
  reader += fscanf(infile, "distance between configs = %zu \n", &delta_config);
  // spatial extend
  reader += fscanf(infile, "Ls = %zu \n", &Ls);
  // temporal extend
  reader += fscanf(infile, "Lt = %zu \n", &Lt);
  // number of eigenvectors
  reader += fscanf(infile, "nb_ev = %zu \n", &nb_ev);
  // starting number of randomvectors
  // TODO: must be changed to allow for several random vectors
  reader += fscanf(infile, "nb_rnd = %zu \n", &nb_rnd);
  // seed for randomvectors
  reader += fscanf(infile, "seed = %i\n", &( seed));
  // verbosity
  reader += fscanf(infile, "verbose = %zu\n", &( verbose));
  // quarktype
  reader += fscanf(infile, "quarktype = %255s\n", readin);
  quarktype.assign(readin);

  // SOURCE --------------------------------------------------------------------
  // type and number of inversions for the source in time ----------------------
  reader += fscanf(infile, "inversion_source_type_t = %255s\n", readin);
  dilution_type_so[0].assign(readin);
  reader += fscanf(infile, "inversion_source_number_t = %zu\n", 
                   &(dilution_size_so[0]));
  // type and number of inversions for the source in eigenvector space ---------
  reader += fscanf(infile, "inversion_source_type_v = %255s\n", readin);
  dilution_type_so[1].assign(readin);
  reader += fscanf(infile, "inversion_source_number_v = %zu\n", 
                   &(dilution_size_so[1]));
  // type and number of inversions for the soure in Dirac space ----------------
  reader += fscanf(infile, "inversion_source_type_d = %255s\n", readin);
  dilution_type_so[2].assign(readin);
  reader += fscanf(infile, "inversion_source_number_d = %zu\n", 
                   &(dilution_size_so[2]));

  // SINK ----------------------------------------------------------------------
  // type and number of dilution vectors for the sink in time ------------------
  reader += fscanf(infile, "inversion_sink_type_t = %255s\n", readin);
  dilution_type_si[0].assign(readin);
  reader += fscanf(infile, "inversion_sink_number_t = %zu\n", 
                   &(dilution_size_si[0]));
  // type and number of dilution vectors for the sink in space -----------------
  reader += fscanf(infile, "inversion_sink_type_s = %255s\n", readin);
  dilution_type_si[1].assign(readin);
  reader += fscanf(infile, "inversion_sink_number_s = %zu\n", 
                   &(dilution_size_si[1]));
  // type and number of dilution vectors for the sink in Dirac space -----------
  reader += fscanf(infile, "inversion_sink_type_d = %255s\n", readin);
  dilution_type_si[2].assign(readin);
  reader += fscanf(infile, "inversion_sink_number_d = %zu\n", 
                   &(dilution_size_si[2]));
  // type and number of dilution vectors for the sink in colour space ----------
  reader += fscanf(infile, "inversion_sink_type_c = %255s\n", readin);
  dilution_type_si[3].assign(readin);
  reader += fscanf(infile, "inversion_sink_number_c = %zu\n", 
                   &(dilution_size_si[3]));

  // output path for the perambulators and the randomvectors
  reader += fscanf(infile, "outpath = %255s\n", readin);
  outpath.assign(readin);
  // input path for the eigensystems
  reader += fscanf(infile, "inpath_ev = %255s\n", readin);
  inpath_ev.assign(readin);

  // close input file
  fclose(infile);
  
  // checking input parameters
  check_input_parameters();
  check_and_create_filenames();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::print_options() {

  std::cout << "config = " <<  config << std::endl;
  std::cout << "Ls = " <<  Ls << ", Lt = " <<  Lt << std::endl;
  std::cout << "nb_ev = " <<  nb_ev << ", nb_rnd = " <<  nb_rnd << std::endl;
  std::cout << "seed = " <<  seed << std::endl;
  std::cout << "verbose = " <<  verbose << std::endl;
  std::cout << "quarktype = " <<  quarktype << std::endl;
  std::cout << "inversion source time: " <<  dilution_type_so[0] << " " 
            <<  dilution_size_so[0] << std::endl;
  std::cout << "inversion source eigenspace: " <<  dilution_type_so[1] << " " 
            <<  dilution_size_so[1] << std::endl;
  std::cout << "inversion source Dirac: " <<  dilution_type_so[2] << " " 
            <<  dilution_size_so[2] << std::endl;
  std::cout << "inversion sink time: " <<  dilution_type_si[0] << " " 
            <<  dilution_size_si[0] << std::endl;
  std::cout << "inversion sink spatial space: " <<  dilution_type_si[1] << " " 
            <<  dilution_size_si[1] << std::endl;
  std::cout << "inversion sink Dirac: " <<  dilution_type_si[2] << " " 
            <<  dilution_size_si[2] << std::endl;
  std::cout << "inversion sink color: " <<  dilution_type_si[3] << " " 
            <<  dilution_size_si[3] << std::endl;
  std::cout << "output path: " <<  outpath << std::endl;
  std::cout << "input path ev: " <<  inpath_ev << std::endl;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------



