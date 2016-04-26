#include "tmLQCD.h"

#include "distillery.h"
#include "input_parms.h"

#include "mpi.h"
#include "omp.h"

extern "C" 
{
#include "global.h"
#include "start.h"
#include "gettime.h"
#include "quda_interface.h"
#include "su3.h"
#include "linalg/square_norm.h"
}

int main(int argc, char *argv[]){

  // MPI initialisation stuff
  // MPI_Init(&argc, &argv);
  int mpi_provided;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &mpi_provided );
  printf("MPI_THREAD_PROVIDED: %d\n", mpi_provided);

  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Group world_group;
  MPI_Comm_group( MPI_COMM_WORLD, &world_group );
  MPI_Comm mpi_comm_world_2;
  MPI_Comm_create( MPI_COMM_WORLD, world_group, &mpi_comm_world_2 );
  

  Eigen::initParallel();
  Eigen::setNbThreads(1);
  
  if(numprocs != 4){
    if(myid == 0){
      std::cout << "Number of processes is: " << numprocs << std::endl;
      std::cout << "Bad choice of MPI processes. Use 4!" << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  LapH::input_parameter param;
  param.parse_input_file(argc, argv);
  // check that everything can be distributed to mpi processes
  if(param.dilution_size_so[2] % numprocs != 0){
    if(myid == 0)
      std::cout << "number of processes must match dilution in Dirac space!\n"
                << std::endl;
    MPI_Finalize();
    return 0;
  } 
  
  int length = 3*4*param.Lt*param.Ls*param.Ls*param.Ls;
  
  // initialisation of the twisted mass stuff
  // set to 1 to make tmLQCD more verbose
  int verbose = 0;
  size_t config = param.config;
  // special version of tmLQCD for QUDA which allocates card no. "myid"
  tmLQCD_invert_init(argc, argv, verbose, myid);
  tmLQCD_read_gauge(config);
  
  // this is for a multithreaded implementation which should be more efficient, essentially hiding
  // part or all of the time spent assembling the perambulator 
  std::complex<double>* source_loc = new std::complex<double>[length];
  std::complex<double>* propagator_loc = new std::complex<double>[length];
  // have QUDA touch memory first
  random_spinor_field_lexic((spinor*)source_loc,0,RN_GAUSS);
  invert_quda_direct((double*)propagator_loc, (double*)source_loc, 1, 0);
   
  std::complex<double>* propagator = new std::complex<double>[length*param.dilution_size_so[2]/numprocs];
  std::complex<double>* source =  new std::complex<double>[length*param.dilution_size_so[2]/numprocs];
  
  // initialisation of distillery
  double atime = gettime();
  LapH::distillery dis;
  dis.initialise(param);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0) printf("Time for distillery initialisation and reading of eigenvectors %.4e seconds\n",gettime()-atime);
  
  if(myid == 0)
    std::cout << "processing config: " << config << "\n" << std::endl;
  // resetting everything
  if((config != param.config)) {
    param.config = config;
    param.nb_rnd = 0;
    if(myid == 0)
      dis.reset_all(param);
  }
  // loop over random vectors
  size_t number_of_random_vectors = 1; // two different random vectors
  for(size_t nb_rnd = 0; nb_rnd < number_of_random_vectors; ++nb_rnd) {
    // loop over all inversions
    omp_set_num_threads(param.peram_gen_omp_num_threads);
    #pragma omp parallel
    {
      // with new CUDA versions, it seems that device memory is dependent on the allocating thread
      // as a result, we need to ensure that only one thread launches the inversion while another
      // takes care of building the perambulator
      int num_threads = omp_get_num_threads();
      int thread_id = omp_get_thread_num();
      // loop over all inversions
      for(size_t dil_t = 0; dil_t < param.dilution_size_so[0]; ++dil_t){
        for(size_t dil_e = 0; dil_e < param.dilution_size_so[1]; ++dil_e){
          if( num_threads == 1 || thread_id == 0 ){
            if(myid == 0) std::cout << "\t\nDoing inversions at: " << dil_t << "\t" << dil_e << "\n" << std::endl;
         
            dis.create_source(dil_t, dil_e, source); 
            MPI_Barrier(mpi_comm_world_2);

            // reordering the sources over processes
            for(int id_source = 0; id_source < numprocs; id_source++){
                MPI_Scatter((double*) (source), 2*length/numprocs, MPI_DOUBLE,
                            (double*) (source_loc+id_source*length/numprocs), 
                            2*length/numprocs, MPI_DOUBLE, id_source, mpi_comm_world_2);
            } 
            MPI_Barrier(mpi_comm_world_2);
          
            invert_quda_direct((double *) propagator_loc, (double *) source_loc, 0, 1);
            MPI_Barrier(mpi_comm_world_2);
          }
          // this barrier is required because otherwise propagator_loc would be accessed from two threads
          // simulteneously
          #pragma omp barrier
          
          if( num_threads == 1 || thread_id == 1 ){
            // reordering the propagators over processes
            for(int id_source = 0; id_source < numprocs; id_source++){
                MPI_Scatter((double*) (propagator_loc), 2*length/numprocs, MPI_DOUBLE,
                            (double*) (propagator+id_source*length/numprocs), 
                            2*length/numprocs, MPI_DOUBLE, id_source, mpi_comm_world_2);
            } 
            MPI_Barrier(mpi_comm_world_2);
          }
          // this barrier is required because otherwise propagator could be accessed from two threads simulateneously
          # pragma omp barrier
          
          // one of the two threads enters this, the other immediately moves to the next iteration
          if( num_threads == 1 || thread_id == 1 ){
            // constructing the perambulator
            for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; dil_d++)
              dis.add_to_perambulator(dil_t, dil_e, dil_d, 
                                      &(propagator[dil_d*length/numprocs]));
          }
          // The MPI fun ends HERE! ;) +++++++++++++++++++++++++++++++++++++++++
        }
      } // end of loop over inversions
    } // OpenMP parallel closing brace

    // creating new random vector and writing perambulator to disk -------------
    dis.write_perambulator_to_disk(); // -------------------------------------
    if(nb_rnd < number_of_random_vectors-1)
      dis.reset_perambulator_and_randomvector(param.seed^(nb_rnd+2), 
                                              nb_rnd+1);
    MPI_Barrier(MPI_COMM_WORLD);
  } // end of loop over random vectors

  tmLQCD_finalise();
  dis.clean();

  MPI_Finalize();
	return 0;
}

