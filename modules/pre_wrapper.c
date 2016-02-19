
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#define MAIN_PROGRAM
#include "config.h"
#include "global.h"
#include "git_hash.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#ifdef MPI
#include "xchange/xchange.h"
#endif
#include <io/utils.h>
#include <io/gauge.h>
#include "read_input.h"
#include "mpi_init.h"
#include "init/init.h"
#include "sighandler.h"
#include "boundary.h"
#include "invert_eo.h"
#include "operator.h"
#include "linalg/convert_eo_to_lexic.h"
#include "pre_wrapper.h"

#ifdef HAVE_GPU
extern void init_mixedsolve_eo(su3** gf);
extern void init_mixedsolve(su3** gf);
extern void finalize_mixedsolve();
extern void init_gpu_fields(int need_momenta);
extern void finalize_gpu_fields();
#include "GPU/cudadefs.h"
  #ifdef TEMPORALGAUGE
     #include "temporalgauge.h" 
   #endif
#endif


void pre_wrapper(int argc, char *argv[]) {

  DUM_DERI = 4; // war 8
  DUM_MATRIX = DUM_DERI + 3; // war 5
  NO_OF_SPINORFIELDS = DUM_MATRIX + 3;

  // in read_input.h
  verbose = 0;
  g_use_clover_flag = 0;

#ifdef MPI
#  ifdef OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif

  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif

  /* Read the input file */
  if( (read_input("invert.input")) != 0) {
    fprintf(stderr, "Could not find input file: invert.input\nAborting...");
  }

#ifdef OMP
  init_openmp();
#endif

  tmlqcd_mpi_init(argc, argv);
  g_dbw2rand = 0;
  for(int j = 0; j < no_operators; j++) if(!operator_list[j].even_odd_flag) even_odd_flag = 0;

#ifndef MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  int j = init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  int j = init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(-1);
  }
  if (even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }
  // define the geometry
  geometry();

  // initialise the operators
  init_operators();

  #ifdef HAVE_GPU
   if(usegpu_flag){
    if(even_odd_flag){
       init_mixedsolve_eo(g_gauge_field);
    }
    else{
      init_mixedsolve(g_gauge_field);
    }
    #ifdef GPU_DOUBLE
      /*init double fields w/o momenta*/
      init_gpu_fields(0);
    #endif
    #ifdef TEMPORALGAUGE
      int retval;
      if((retval=init_temporalgauge(VOLUME, g_gauge_field)) !=0){
	      if(g_proc_id == 0) printf("Error while initializing temporal gauge. Aborting...\n");   
	      exit(200);
      }
    #endif
   }//usegpu_flag
  #endif  

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }
  if (g_sloppy_precision_flag == 1) {
    j = init_dirac_halfspinor32();
    if (j != 0) {
      fprintf(stderr, "Not enough memory for 32-bit halffield! Aborting...\n");
      exit(-1);
    }
  }
#  if (defined _PERSISTENT)
  if (even_odd_flag)
    init_xchange_halffield();
#  endif
#endif

  return;
}

int tmLQCD_read_gauge(const int nconfig) {
  char conf_filename[500];
  sprintf(conf_filename, "%s.%.4d", gauge_input_filename, nconfig);
  int j=0;
  if (g_cart_id == 0) {
    printf("#\n# Trying to read gauge field from file %s.\n",
	   conf_filename);
    fflush(stdout);
  }
  if( (j = read_gauge_field(conf_filename)) !=0) {
//  if( (j = read_gauge_field(conf_filename,g_gauge_field)) !=0) {
    fprintf(stderr, "Error %d while reading gauge field from %s\n Aborting...\n", j, conf_filename);
    exit(-2);
  }
  if (g_cart_id == 0) {
    printf("# Finished reading gauge field.\n");
    fflush(stdout);
  }
#ifdef MPI
  xchange_gauge(g_gauge_field);
#endif
  return(0);
}


int invert(double * propagator, double * source) {
  unsigned int op_id = 0;
  unsigned int write_prop = 0;
  unsigned int index_start = 0;
  g_mu = 0.;

  operator_list[op_id].sr0 = g_spinor_field[0];
  operator_list[op_id].sr1 = g_spinor_field[1];
  operator_list[op_id].prop0 = g_spinor_field[2];
  operator_list[op_id].prop1 = g_spinor_field[3];

  // convert to even/odd geometry
  convert_lexic_to_eo(operator_list[op_id].sr0, operator_list[op_id].sr1, (spinor*) source);
  
  operator_list[op_id].inverter(op_id, index_start, write_prop);
  
  // convert back to lexicographic geometry
  convert_eo_to_lexic((spinor*) propagator, operator_list[op_id].prop0, operator_list[op_id].prop1);


  return(0);
}
