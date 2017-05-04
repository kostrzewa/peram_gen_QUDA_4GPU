#include "./distillery.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::initialise(const LapH::input_parameter& in_param) {

  param = in_param;
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;

  int numprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // Initializing memory for eigenvectors, perambulator and random vector
  V = new Eigen::MatrixXcd[Lt/numprocs];
  perambulator = new Eigen::MatrixXcd[Lt/numprocs];
  for(size_t t = 0; t < Lt/numprocs; ++t){
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
    perambulator[t] = Eigen::MatrixXcd::Zero( 4 * param.nb_ev,
                                                  param.dilution_size_so[0] * 
                                                  param.dilution_size_so[1] * 
                                                  param.dilution_size_so[2] );
  }
  random_vector = new Eigen::VectorXcd[Lt];
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);

  // reading eigenvectors from disk
  read_eigenvectors(); 
  // generating random vector
  set_random_vector();
  // is everything allocated?
  memory_allocated = true;
  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::clean() {

  int numprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // matrices are resized to size 0
  for(size_t t = 0; t < param.Lt/numprocs; ++t){
    (V[t]).resize(0, 0);
    (perambulator[t]).resize(0, 0);
    (random_vector[t]).resize(0);
  }
  for(size_t t = 0; t < param.Lt; ++t){
    (random_vector[t]).resize(0);
  }
  delete[] V;
  delete[] perambulator;
  delete[] random_vector;
  V = NULL;
  perambulator = NULL;
  random_vector = NULL;

  memory_allocated = false;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_all(const LapH::input_parameter& in_param){

  param = in_param;
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;

  int numprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  for(size_t t = 0; t < Lt/numprocs; ++t){
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
    perambulator[t] = Eigen::MatrixXcd::Zero( 4 * param.nb_ev,
                                                  param.dilution_size_so[0] * 
                                                  param.dilution_size_so[1] * 
                                                  param.dilution_size_so[2] );
  }
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);

  read_eigenvectors();  
  set_random_vector();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_perambulator_and_randomvector(const int seed,
                                            const size_t random_vector_number){

  int numprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // setting perambulator and random vector to zero 
  for(size_t t = 0; t < param.Lt/numprocs; ++t)
    perambulator[t] = Eigen::MatrixXcd::Zero( 4 * param.nb_ev,
                                                  param.dilution_size_so[0] * 
                                                  param.dilution_size_so[1] * 
                                                  param.dilution_size_so[2] );
  for(size_t t = 0; t < param.Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);

  // resetting random seed and generating new random vector
  param.nb_rnd = random_vector_number;
  param.seed = seed; 
  set_random_vector();

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Creates the dilution lookup table. 
// input:  dilution_size  -> size of dilution, e.g. block size or interlace size
//         dilution_entry -> index of the dilution
//         type           -> just the type of the dilution
// output: lookup         -> the final lookup table for dilution
static void create_dilution_lookup(const size_t nb_of_nonzero_entries, 
                                   const size_t nb_inversions, 
                                   const size_t dilution_entry, 
                                   const std::string type, size_t lookup[]){
  if(type == "B" || type == "N")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = nb_of_nonzero_entries * dilution_entry + i;
  else if(type == "I" || type == "F")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = dilution_entry + i * nb_inversions;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::create_source(const size_t dil_t, const size_t dil_e,
                                     std::complex<double>* source) {
  
  if(dil_t >= param.dilution_size_so[0] ||
     dil_e >= param.dilution_size_so[1] ){
        std::cout << "dilution is out of bounds in \"create_source\"" <<
        std::endl;
        exit(0);
  }

  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t Vs = Ls*Ls*Ls;
  const size_t dim_row = Vs*3;

  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // setting source to zero 
  for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d)
    for(size_t i = 0; i < 12*Vs*Lt/numprocs; ++i)
      source[dil_d*12*Vs*Lt/numprocs + i] = {0.0, 0.0};

  // indices of timeslices with non-zero entries
  size_t nb_of_nonzero_t = Lt/param.dilution_size_so[0];  
  // TODO: think about smart pointer here!
  size_t* t_index = new size_t[nb_of_nonzero_t];
  create_dilution_lookup(nb_of_nonzero_t, param.dilution_size_so[0], 
                         dil_t, param.dilution_type_so[0], t_index);

  // indices of eigenvectors to be combined
  size_t nb_of_ev_combined = number_of_eigen_vec/param.dilution_size_so[1];
  size_t* ev_index = new size_t[nb_of_ev_combined];
  create_dilution_lookup(nb_of_ev_combined, param.dilution_size_so[1], 
                         dil_e, param.dilution_type_so[1], ev_index);

  // indices of Dirac components to be combined 
  // TODO: This is needed only once - could be part of class
  size_t nb_of_dirac_combined = 4/param.dilution_size_so[2];
  size_t** d_index = new size_t*[param.dilution_size_so[2]];
  for(size_t i = 0; i < param.dilution_size_so[2]; ++i)
    d_index[i] = new size_t[nb_of_dirac_combined];
  for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d)
    create_dilution_lookup(nb_of_dirac_combined, param.dilution_size_so[2],
                           dil_d, param.dilution_type_so[2], d_index[dil_d]);

  // creating the source 
  // running over nonzero timeslices
  for(size_t t = 0; t < nb_of_nonzero_t; ++t){ 
    // intermidiate memory
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(dim_row, 4);
    size_t time = t_index[t];                 // helper index
    if((int) (time / (Lt/numprocs)) == myid) {
      size_t time_id = time % (Lt/numprocs);  // helper index
      // building source on one timeslice
      for(size_t ev = 0; ev < nb_of_ev_combined; ++ev){
        size_t ev_h = ev_index[ev] * 4; // helper index 
        for(size_t d = 0; d < 4; ++d){
          S.col(d) += random_vector[time](ev_h+d) *
                          (V[time_id]).col(ev_index[ev]); 
        }
      }
      // copy the created source into output array
      size_t t_h = time_id*Vs; // helper index
      for(size_t x = 0; x < Vs; ++x){
        size_t x_h  = (t_h + x)*12;    // helper index
        size_t x_h2 = x*3;            // helper index
        for(size_t d2 = 0; d2 < param.dilution_size_so[2]; ++d2){
          for(size_t d3 = 0; d3 < nb_of_dirac_combined; ++d3){
            size_t d_h = x_h + d_index[d2][d3]*3; // helper index
            for(size_t c = 0; c < 3; ++c){
              source[d2*12*Vs*Lt/numprocs + d_h + c] = 
                                          S(x_h2 + c, d_index[d2][d3]);
            }
          }
        }    
      }
    }
  } // end of loop over nonzero timeslices
  MPI_Barrier(MPI_COMM_WORLD);
  // freeing memory
  delete[] t_index;
  delete[] ev_index;
  for(size_t i = 0; i < param.dilution_size_so[2]; ++i)
    delete[] d_index[i];
  delete[] d_index;
  t_index = NULL;
  ev_index = NULL;
  d_index = NULL;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::add_to_perambulator(const size_t dil_t, 
                          const size_t dil_e, const size_t dil_d,
                          const std::complex<double>* const propagator) { 
 
  if(dil_t >= param.dilution_size_so[0] ||
     dil_e >= param.dilution_size_so[1] ||
     dil_d >= param.dilution_size_so[2] ){
        std::cout << "dilution is out of bounds in \"add_to_perambulator\"" <<
        std::endl;
        exit(0);
  }
 
  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t vec_length = Ls*Ls*Ls*3;
  const size_t nb_ev = param.nb_ev;
  // computing the column index of perambulator
  const size_t col = dil_t*param.dilution_size_so[1]*param.dilution_size_so[2] +
                     dil_e*param.dilution_size_so[2] + dil_d;

  int numprocs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // intermediate memory for multiplication with V^dagger
  Eigen::MatrixXcd vec = Eigen::MatrixXcd::Zero(vec_length, 4);

  // running over sink time index
  for(size_t t = 0; t < Lt/numprocs; ++t){
    // running over all indices on one timeslice 
    // -> resorting propagator and copyying it into vec
    size_t t_h = 12*Ls*Ls*Ls*t; // helper index
    for(size_t x = 0; x < Ls*Ls*Ls; ++x){ // spatial
      size_t x_h = t_h + 12*x;  // helper index
      for(size_t d = 0; d < 4; ++d){      // Dirac
        size_t d_h = x_h + 3*d; // helper index
        for(size_t c = 0; c < 3; ++c){    // colour
           vec(3*x + c, d) = propagator[d_h + c];
        }
      }
    }
  
    // multiplication with V^dagger and storing result in perambulator
    Eigen::MatrixXcd vec1 = Eigen::MatrixXcd::Zero(nb_ev, 4);
    for(size_t d = 0; d < 4; ++d)    
      vec1 = V[t].adjoint()*vec;
    for(size_t c = 0; c < nb_ev; ++c)  
      for(size_t d = 0; d < 4; ++d)    
        (perambulator[t])(4*c+d, col) = vec1(c, d);

  } // end of loop over sink time index

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::write_perambulator_to_disk() {

  char outfile[400];
  FILE *fp = NULL;
  const size_t Lt = param.Lt;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;
  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  const size_t number_of_inversions = param.dilution_size_so[0] *
                                      param.dilution_size_so[1] *
                                      param.dilution_size_so[2] ;
  const size_t size_perambulator_entry = number_of_inversions * Lt * 4 * 
                                         number_of_eigen_vec / numprocs;

  // memory for reading perambulators
  std::complex<double>* perambulator_write =
      new std::complex<double>[size_perambulator_entry];
  // copy perambulator into writing array
  for(size_t t = 0; t < Lt/numprocs; ++t){
    size_t t_h = t*4*number_of_eigen_vec; // helper index
    for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i){
      size_t row_i_h = (row_i + t_h) * number_of_inversions; // helper index
      for(size_t col_i = 0; col_i < number_of_inversions; ++col_i){
        perambulator_write[row_i_h + col_i] = perambulator[t](row_i, col_i);
      }
    }
  }

  // data path
  std::string filename = param.outpath + "/perambulator";
  if(verbose) printf("writing perambulators from files:\n");
  else printf("\twriting perambulator\n");
  // create perambulator file name
  sprintf(outfile,
      "%s.rndvecnb%02d.%s.Tso%s%04d.Vso%s%04d.Dso%s%01d.Tsi%s%04d.Ssi%s%04d" \
      ".Dsi%s%01d.Csi%s%01d.smeared%01d.%05d",
      filename.c_str(), (int) param.nb_rnd, param.quarktype.c_str(), 
      param.dilution_type_so[0].c_str(), (int) param.dilution_size_so[0], 
      param.dilution_type_so[1].c_str(), (int) param.dilution_size_so[1], 
      param.dilution_type_so[2].c_str(), (int) param.dilution_size_so[2],
      param.dilution_type_si[0].c_str(), (int) param.dilution_size_si[0], 
      param.dilution_type_si[1].c_str(), (int) param.dilution_size_si[1], 
      param.dilution_type_si[2].c_str(), (int) param.dilution_size_si[2], 
      param.dilution_type_si[3].c_str(), (int) param.dilution_size_si[3], 
      0, (int) param.config); // we don't have support for sink smearing in this version
  // writing data
  for(int id = 0; id < numprocs; id++){
    if(myid == id){

      if(myid == 0){ // is needed if file already exists
        if((fp = fopen(outfile, "wb")) == NULL){
          std::cout << "failed to open file: " << outfile << "\n" << std::endl;
          MPI_Finalize();
          exit(0);
        }
      }
      else{
        if((fp = fopen(outfile, "ab")) == NULL){
          std::cout << "failed to open file: " << outfile << "\n" << std::endl;
          MPI_Finalize();
          exit(0);
        }
      }
      if(verbose) printf("\twrite file: %s\n", outfile);
      fwrite(perambulator_write, sizeof(std::complex<double>),
          size_perambulator_entry, fp);
      fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);  
  }

  delete[] perambulator_write;

}
////////////////////////////////////////////////////////////////////////////////
/////////////////////// private methods ////////////////////////////////////////
void LapH::distillery::set_random_vector() {

  double sqrt2 = 0.5*sqrt(2.0);
  double re, im;
  rlxs_init(2, param.seed^param.config);
  int rnd_length = 2*param.Lt*param.nb_ev*4;
  float* rnd = new float[rnd_length];
  ranlxs(rnd, rnd_length);

  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // generating a Z_2 source
  for(size_t t = 0; t < param.Lt; ++t ){ 
    random_vector[t].setRandom(4 * param.nb_ev);
    for(size_t i = 0; i < 4 * param.nb_ev; ++i){
      if (rnd[8 * param.nb_ev * t + 2*i] < 0.5)
        re = -sqrt2;
      else 
        re = sqrt2;
      if (rnd[8 * param.nb_ev * t + 2*i + 1] < 0.5)
        im = -sqrt2;
      else 
        im = sqrt2;
      random_vector[t](i) = {re, im};
    }
  }


  // writing random vector to disc
  if(myid == 0)
    write_random_vector_to_disk();

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::read_random_vector() {

  char infile[400];
  FILE *fp = NULL;
  const int Lt = param.Lt;
  const int verbose = param.verbose;
  const int number_of_rnd_vec = 1;
  const int number_of_eigen_vec = param.nb_ev;
  const int rnd_vec_length = number_of_eigen_vec * 4;
  // memory for reading random vectors
  std::complex<double>* rnd_vec_read =
      new std::complex<double>[Lt*rnd_vec_length];
  std::string filename = "/data/LapHs/peram_gen/new/test/randomvector";

  if(verbose) printf("reading random vectors from files:\n");
  else printf("\treading random vectors:");
  int check_read_in = 0;

  for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
    // read random vector
    sprintf(infile, "%s.%03d.u.Tf.%04d", filename.c_str(), 0, 1000);
    if(verbose) printf("\tread file: %s\n", infile);
    if((fp = fopen(infile, "rb")) == NULL){
      std::cout << "failed to open file: " << infile << "\n" << std::endl;
      exit(0);
    }   
    check_read_in += fread(rnd_vec_read, sizeof(std::complex<double>),
        Lt*rnd_vec_length, fp);
    // copy into matrix structure
    for(int t = 0; t < Lt; t++){
      for(int row_i = 0; row_i < rnd_vec_length; ++row_i){
        random_vector[t](row_i) = rnd_vec_read[t * rnd_vec_length + row_i];
      }   
    }
    fclose(fp);
  }   
  delete[] rnd_vec_read;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::write_random_vector_to_disk(){

  char outfile[400];
  FILE *fp = NULL;
  const size_t Lt = param.Lt;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t rnd_vec_length = Lt * number_of_eigen_vec * 4;

  // copy from eigen structure into intermediate memory
  std::complex<double>* rnd_vec_write =
      new std::complex<double>[rnd_vec_length];
  for(size_t t = 0; t < Lt; ++t)
    for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i)
      rnd_vec_write[row_i + t * rnd_vec_length/Lt] = random_vector[t](row_i);

  // creating name and path of outfile
  std::string filename = param.outpath + "/randomvector";
  if(verbose) printf("writing random vector to files:\n");
  else printf("\twriting random vector\n");
  int check_read_in = 0;
  // TODO: Must be changed to write stochastic sink random vectors
  sprintf(outfile, "%s.rndvecnb%02d.%s.nbev%04d.%04d", filename.c_str(), 
                    (int) param.nb_rnd, param.quarktype.c_str(), 
                    (int) param.nb_ev, (int) param.config);
  if(verbose) printf("\twrite to file: %s\n", outfile);
  if((fp = fopen(outfile, "wb")) == NULL){
    std::cout << "failed to open file to write random vector: " 
              << outfile << "\n" << std::endl;
    exit(0);
  }   
  check_read_in += fwrite(rnd_vec_write, sizeof(std::complex<double>),
                          rnd_vec_length, fp);
  if(check_read_in != (int) rnd_vec_length){
    std::cout << "failed to write random vector: "
              << outfile << "\n" << std::endl;
    exit(0);
  }   
  fclose(fp);
  // delete intermediate memory
  delete[] rnd_vec_write;

} 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::read_eigenvectors(){

  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;

  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];

  if(verbose) printf("reading eigen vectors from files:\n");
  else printf("\treading eigenvectors\n");
  fflush(stdout);

  size_t start = myid*Lt/numprocs;
  size_t end = start + Lt/numprocs;
  size_t counter = 0;
  for(size_t t = start; t < end; ++t){
    //setting up file
    char name[200];
    sprintf(name, "%s/eigenvectors.%04d.%03d", 
                  param.inpath_ev.c_str(), (int) param.config, (int) t);
    if(verbose) std::cout << "Reading file: " << name << std::endl;
    std::ifstream infile(name, std::ifstream::binary);
    if (infile) {
      for (size_t nev = 0; nev < number_of_eigen_vec; ++nev) {
        infile.read( (char*) eigen_vec, 2*dim_row*sizeof(double));
        for(size_t nrow = 0; nrow < dim_row; ++nrow)
          (V[counter])(nrow, nev) = eigen_vec[nrow];
      }
    }
    else {
      std::cout << "eigenvector file does not exist!!!\n" << std::endl;
      exit(0);
    }
    infile.close();

    // small test of trace and sum over the eigen vector matrix!
    if(verbose){
      std::cout << "trace of V^d*V on t = " << counter << ":\t"
          << (V[counter].adjoint() * V[counter]).trace() << std::endl;
      std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"
          << (V[counter].adjoint() * V[counter]).sum() << std::endl;
    }
    counter++;
  }
  delete[] eigen_vec;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------












