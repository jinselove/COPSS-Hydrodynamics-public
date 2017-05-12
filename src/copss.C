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

	comm_in = init.comm();

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
  if(comm_in.rank()==0)
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
  if(comm_in.rank()==0)
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
  const GetPot tmp(control_file);
  input_file = tmp;
  this -> read_test_name();
  this -> read_physical_parameter();
  this -> read_particle_parameter();
  this -> read_geometry();
} // end read_data function


//====================================================================
void Copss::read_test_name()
  {
    test_name = input_file("test_name", "validation");

    if(comm_in.rank() == 0){
      printf("##########################################################\n"
             "#                       system_name                       \n"
             "##########################################################\n\n"
             "-----------> system_name: %s\n", test_name.c_str() );
    }
}

  /*
   * Read physical parameters 
   */
void Copss::read_physical_parameter()
{
  T = input_file("temperature", 297);// K
  kBT    = kB * T; //(N*um)
  viscosity            = input_file("viscosity", 1.0); // viscosity (cP = N*s/um^2)
  Rb                   = input_file("radius", 0.10); // radius of the bead (um)
  drag_c      = 6.*PI*viscosity*Rb;    // Drag coefficient (N*s/um)
  Db          = kBT/drag_c;     // diffusivity of a bead (um^2/s)  

  tc   = drag_c*Rb*Rb/kBT;       // diffusion time (s)
  uc   = kBT/(drag_c*Rb);        // characteristic velocity (um/s)
  fc   = kBT/Rb;                 // characteristic force (N)
  muc  = 1./(6.*PI);             // non-dimensional viscosity

  // print out physical parameters information
  if(comm_in.rank() == 0){
    printf(" ##########################################################\n"
           " #                  System Physical Parameters             \n"
           " ##########################################################\n\n"
           "   temperature           T   =  %.6e (K)\n"
           "   viscosity             mu  =  %.6e (cP = N*s/um^2)\n"
           "   Energy unit           kBT =  %.6e (N*um = N*um)\n"
           "   Radius of the bead     a  =  %.6e (um)\n"
           "   bead diffusivity      Db  =  %.6e (um^2/s)\n"
           "   HI Drag coefficient  zeta = 6*PI*mu*a =  %.6e (N*s/um)\n"
           "   ksi = sqrt(PI)/(3a)       =  %.6e (1/um)\n"
           "   ------------> The characteristic variables:\n"
           "   characteristic time          = %.6e (s)\n"
           "   characteristic velocity      = %.6e (um/s)\n"
           "   characteristic force         = %.6e (N)\n",
           T, viscosity, kBT, Rb, Db, drag_c, std::sqrt(PI)/(3.*Rb), tc, uc, fc);
  } // end if (comm_in.rank() == 0)
}// end read_physical_parameter()

  /*
   * Read Geometry infomation
   */

void Copss::read_geometry()
{
  dim = input_file("dimension", 3);
  //=============== wall type and wall params
  wall_type = input_file("wall_type", "not_defined");
  wall_params.resize(input_file.vector_variable_size(wall_type));
  if(wall_type != "not_defined"){
    for (unsigned int j = 0; j < wall_params.size(); j++){
      wall_params[j] = input_file(wall_type,0.0,j);
    }
  }
  else{
    error_msg = "wall_type undefined; please check the wall_type definition in control file (1. wall_type; 2. wall_params)";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  //=============== periodicity
  periodicity.resize(input_file.vector_variable_size("periodicity"));
  for (unsigned int i=0; i < periodicity.size(); i++){ periodicity[i] = input_file("periodicity", false, i); }
  if(periodicity[0]==true and periodicity[1] == true and periodicity[2]==true){
    error_msg = "warning: The box cannot be periodic on all directions at the same time. (required by FEM)";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  //============== inlet 
  inlet.resize(input_file.vector_variable_size("inlet"));
  for (unsigned int i=0; i < inlet.size(); i++){ 
     inlet[i] = input_file("inlet", false, i);
     if(inlet[i]==true and periodicity[i]==true) {
      error_msg = "warning: A inlet direction has to be non-periodicity";
      PMToolBox::output_message(error_msg,comm_in);
      libmesh_error();
     }
  }
  //============== inlet pressure
  inlet_pressure.resize(input_file.vector_variable_size("inlet_pressure"));
  for (unsigned int i=0; i < inlet_pressure.size(); i++){ inlet_pressure[i] = input_file("inlet_pressure", 0, i); }

  if(comm_in.rank() == 0){
    printf ("##########################################################\n"
           "#                  Geometry information                   \n" 
           "##########################################################\n\n"
           "  Dimension: %d \n"
           "  Wall type: %s \n"
           "  Wall parameters: ",
           dim, wall_type.c_str());
    for (int i = 0; i < wall_params.size(); ++i){
      printf("%f, ", wall_params[i]);
    }
    printf("\n");

    printf("  Periodicity of the box: ");
    for (int i = 0; i < dim; i++){
      printf("%s, ", periodicity[i] ? "Ture" : "False"); 
    }
    printf("\n");

    printf("  Inlet/Outlet of the box: ");
    for (int i = 0; i < dim; i++){
      printf("%s (pressure = %.4e), ", periodicity[i] ? "Ture" : "False", inlet_pressure[i]); 
    }
    printf("\n");
  } // end if comm_in.rank() == 0
} // end read_geometry()

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




