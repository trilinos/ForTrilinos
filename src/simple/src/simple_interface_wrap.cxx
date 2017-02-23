#include "simple_interface.hpp"

#ifdef __cplusplus
extern "C" {
#endif
// void swigc_init__SWIG_0(void* farg1) {
  // ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();
//
  // MPI_Comm arg1 ;
//
  // arg1 = *(MPI_Comm *)(farg1);
  // si.init(arg1);
// }


void swigc_init__SWIG_1() {
  ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();

  si.init();
}


void swigc_setup_matrix(int* farg1, int* farg2, int* farg3, int* farg4, int* farg5, double* farg6) {
  ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();


  int arg1 ;
  int *arg2 = (int *) 0 ;
  int *arg3 = (int *) 0 ;
  int arg4 ;
  int *arg5 = (int *) 0 ;
  double *arg6 = (double *) 0 ;

  arg1 = *farg1;
  arg2 = farg2;
  arg3 = farg3;
  arg4 = *farg4;
  arg5 = farg5;
  arg6 = farg6;

  // FIXME: do a proper conversion
  for (int i = 0; i < arg1; i++)
    arg2[i]--;
  for (int i = 0; i < arg4; i++)
    arg5[i]--;

  si.setup_matrix(arg1,(int const *)arg2,(int const *)arg3,arg4,(int const *)arg5,(double const *)arg6);
}


// void swigc_setup_solver(void* farg1) {
  // ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();
//
  // Teuchos::RCP< Teuchos::ParameterList > *arg1 = 0 ;
//
  // arg1 = (Teuchos::RCP< Teuchos::ParameterList > *)(farg1);
  // si.setup_solver((Teuchos::RCP< Teuchos::ParameterList > const &)*arg1);
// }


void swigc_solve(int* farg1, double* farg2, double* farg3) {
  ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();

  int arg1 ;
  double *arg2 = (double *) 0 ;
  double *arg3 = (double *) 0 ;

  arg1 = *farg1;
  arg2 = farg2;
  arg3 = farg3;
  si.solve(arg1,(double const *)arg2,arg3);
}


void swigc_finalize() {
  ForTrilinos::SimpleInterface& si = ForTrilinos::SimpleInterface::getInstance();

  si.finalize();
}


#ifdef __cplusplus
}
#endif
