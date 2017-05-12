// Parallel Finite Element-General Geometry Ewald-like Method.
// Copyright (C) 2015-2016 Xujun Zhao, Jiyuan Li, Xikai Jiang

// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.


// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.


// You should have received a copy of the GNU General Public
// License along with this code; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include "copss.h"

namespace libMesh
{

//==========================================================================
Copss::Copss(int argc,char** argv)
{
	LibMeshInit init(argc, argv);

	this -> check_libmesh();

	_comm_in = init.comm();

}


int Copss::check_libmesh(){

  // This example NaNs with the Eigen sparse linear solvers and Trilinos solvers,
  // but should work OK with either PETSc or Laspack.
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS,
                           "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS,
                           "--enable-petsc or --enable-laspack");
  
  // check libmesh enables
  #ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
  #endif
  
  #ifndef LIBMESH_HAVE_SLEPC
  libmesh_example_requires(false, "--enable-slepc");
  #endif
  
  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D/3D support");

  return 0;
}


//==========================================================================
void Copss::start_time(struct tm * timeinfo){
  if(_comm_in.rank()==0)
  {
    printf("\n");
    printf("---------------------------------------------------------------------\n");
    printf("The program starts. \n");
    printf("The current date/time is: %s", asctime (timeinfo) );
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
}


//==========================================================================
void Copss::end_time(struct tm * timeinfo){
  if(_comm_in.rank()==0)
  {
    printf("\n");
    printf("---------------------------------------------------------------------\n");
    printf("The current date/time is: %s", asctime (timeinfo) );
    printf("The program ends. \n");
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }	
}


//====================================================================
  void Copss::read_data(std::string control_file)
  {
    //GetPot tmp(control_file);
    //_test_name = tmp("test_name" , "validation");
    const GetPot tmp(control_file);
    _input_file = tmp;
    this -> read_test_name();
  } // end read_data function


//====================================================================
  void Copss::read_test_name()
  {
    _test_name = _input_file("test_name", "validation");

    if(_comm_in.rank() == 0){
      printf("##########################################################\n"
             "#                       system_name                       \n"
             "##########################################################\n\n"
             "-----------> system_name: %s\n", _test_name.c_str() );
    }
  }

  /*
   * Read physical parameters 
   */
  void Copss::read_physical_parameters(){}

  /*
   * Read Geometry infomation
   */

  void Copss::read_geometry(){}

  /*
   * Read mesh
   */
  void Copss::read_mesh(){}

  /*
   * read force types
   */
  void Copss::read_force(){}

  /*
   * read GGEM info
   */
  void Copss::read_ggem(){}

  /*
   * read Stokes Solver  
   */
  void Copss::read_stokes_solver(){}

  /*
   * read Chebyshev info
   */
  void Copss::read_chebyshev(){}

  /*
   * read run time info 
   */
  void Copss::read_run(){}


}




