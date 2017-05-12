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
  this -> read_test_info();
  this -> read_physical_info();
  this -> read_particle_info();
  this -> read_geometry_info();
  this -> read_mesh_info();
  this -> read_force_info();
  this -> read_ggem_info();
  this -> read_stokes_solver_info();
  this -> read_chebyshev_info();
} // end read_data function


//====================================================================
void Copss::read_test_info()
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
void Copss::read_physical_info()
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

void Copss::read_geometry_info()
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
void Copss::read_mesh_info(){
  if(comm_in.rank() == 0){
  printf("##########################################################\n"
         "#                   Mesh information                      \n"
         "##########################################################\n\n");
  }

  generate_mesh = input_file("generate_mesh", false); 
  if (generate_mesh){
    n_mesh.resize(input_file.vector_variable_size("n_mesh"));
    for (unsigned int i=0; i < n_mesh.size(); i++){ n_mesh[i] = input_file("n_mesh", 10, i); }
    if(comm_in.rank() == 0){
      printf(" Generate Mesh:  n_mesh = ");
      for (int i = 0 ; i < dim; i++) printf("%.2e, ",n_mesh[i]);
      printf("\n");
    } // end if comm_in.rank() == 0
  }
  else{
    domain_mesh_file = input_file("domain_mesh_file" , "nothing");
    if(comm_in.rank() == 0){
      printf(" Load mesh file from = %s\n", domain_mesh_file.c_str());
    } // end if comm_in.rank() == 0  
  } // end else

} // end read_mesh_info()


  /*
   * read force types
   */
void Copss::read_force_info(){

  // read particle-particle force types
  num_pp_forces = input_file.vector_variable_size("particle_particle_force_types");
  pp_force_types.resize(num_pp_forces);
  pp_forces.resize(num_pp_forces);
  for (unsigned int i=0; i < num_pp_forces; i++){    
    pp_force_types[i] = input_file("particle_particle_force_types", "nothing", i);
    std::vector<Real> params(input_file.vector_variable_size(pp_force_types[i]));
    if(pp_force_types[i] != "nothing"){
      for (unsigned int j = 0; j < params.size(); j++){
        params[j] = input_file(pp_force_types[i],0.0,j);
      }
    }
    pp_forces[i].first = pp_force_types[i];
    pp_forces[i].second = params;
  }

    // read particle-wall force types
  num_pw_forces = input_file.vector_variable_size("particle_wall_force_types");
  pw_force_types.resize(num_pw_forces);
  pw_forces.resize(num_pw_forces);
  for (unsigned int i=0; i < num_pw_forces; i++){    
    pw_force_types[i] = input_file("particle_wall_force_types", "nothing" , i);
    std::vector<Real> params(input_file.vector_variable_size(pw_force_types[i]));
    if(pw_force_types[i] != "nothing"){
      for (unsigned int j = 0; j < params.size(); j++){
        params[j] = input_file(pw_force_types[i],0.0,j);
      }
    }
    pw_forces[i].first = pw_force_types[i];
    pw_forces[i].second = params;
  } 
  if(comm_in.rank() == 0){
    printf("##########################################################\n"
          "#    Force information (particle-particle)                \n"
          "##########################################################\n\n");
    for (int i = 0; i < num_pp_forces; i++){
      printf ("  ");
      printf("%s  ", pp_forces[i].first.c_str());
      for (int j = 0; j < pp_forces[i].second.size(); j++){
        printf("%.6e  ", pp_forces[i].second[j]);      
      }
      printf("\n");
    }
    printf("##########################################################\n"
          "#    Force information (particle-wall)                \n"
          "##########################################################\n\n");
    for (int i = 0; i < num_pw_forces; i++){
      printf ("  ");
      printf("%s  ", pw_forces[i].first.c_str());
      for (int j = 0; j < pw_forces[i].second.size(); j++){
        printf("%.6e  ", pw_forces[i].second[j]);      
      }
      printf("\n");
    }
  } // end if comm_in.rank() == 0
} // end read_force_info()

/*
 * read GGEM info
 */
void Copss::read_ggem_info(){
  alpha                = input_file("alpha", 0.1);
  if(comm_in.rank() == 0){
    printf("##########################################################\n"
           "#                 GGEM information                       \n"
           "##########################################################\n\n"
           "  The smoothing parameter in GGEM alpha = %.4e\n" 
           "  Recommend meshsize <= %.4e\n", alpha, 1./(std::sqrt(2)*alpha));
  }
}

/*
 * read Stokes Solver  
 */
void Copss::read_stokes_solver_info(){
  max_linear_iterations = input_file("max_linear_iterations", 100);
  linear_solver_rtol   = input_file("linear_solver_rtol", 1E-6);
  linear_solver_atol   = input_file("linear_solver_atol", 1E-6);
  user_defined_pc            = input_file("user_defined_pc", true);
  schur_user_ksp       = input_file("schur_user_ksp", false);
  schur_user_ksp_rtol  = input_file("schur_user_ksp_rtol", 1E-6);
  schur_user_ksp_atol  = input_file("schur_user_ksp_atol", 1E-6);
  schur_pc_type = input_file("schur_pc_type", "SMp");
  stokes_solver_type = input_file("stokes_solver", "superLU_dist");
  if(stokes_solver_type=="superLU_dist") {
    solver_type = superLU_dist;
    user_defined_pc = false;
  }
  else if(stokes_solver_type=="field_split") {
    solver_type = field_split;
    user_defined_pc = true;
  }
  else {
    solver_type = user_define;
  }  

  if(comm_in.rank() == 0){
    printf("##########################################################\n"
           "#                 Solver information                     \n"
           "##########################################################\n\n"
           "  Stokes solver type = %s\n", stokes_solver_type.c_str());
    if (stokes_solver_type=="field_split"){
      printf("  FieldSplit Schur Complement Reduction Solver \n"
             "  schur_pc_type = %s\n", schur_pc_type.c_str());
        if(schur_user_ksp){
              printf(" user defined KSP is used for Schur Complement!\n"
                     " KSP rel tolerance for Schur Complement solver is = %.4e\n",
                     " KSP abs tolerance for Schur Complement solver is = %.4e\n",
              schur_user_ksp_rtol,schur_user_ksp_atol);
          }// end if(schur_user_ksp)
    }// end if (stokes_solver_type == "field_split")
  } // end if (comm_in.rank() == 0)
}// end read_stokes_solver_info()

/*
 * read Chebyshev info
 */
void Copss::read_chebyshev_info(){
    max_n_cheb = input_file("max_n_cheb", 10);
    tol_cheb = input_file("tol_cheb", 0.1);
    eig_factor = input_file("eig_factor", 1.05);
    tol_eigen = input_file("tol_eigen", 0.01);

    // Initially set compute_eigen flag to be true
    compute_eigen = true;

    // print out information
    if(comm_in.rank() == 0){
      printf("##########################################################\n"
             "#   Chebyshev information (only needed by brownian System)   \n"
             "##########################################################\n\n"  
             "  Initially compute_eigen flag is set to be True !!!!\n"
             "  max number of chebyshev polynomial = %d\n"
             "  tolerance of chebyshev polynomial = %.4e\n"
             "  factor of eigenvalues range = %.4e\n"
             "  tolerance of eigenvalues convergence = %.4e\n",  
             max_n_cheb, tol_cheb, eig_factor, tol_eigen);
    } // end if (comm_in.rank() == 0)

} // end read_chebyshev_info()


/*
 * read run time info 
 */
void Copss::read_run_info(){}


}




