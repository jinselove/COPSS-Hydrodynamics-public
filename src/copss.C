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
Copss::Copss(const CopssInit& init)
{
  //mpi_initialized = 1;

  comm_in = init.comm();

	this -> check_libmesh();

}

Copss::~Copss()
{

  delete mesh;
  delete pm_periodic_boundary;
  delete force_field;
  delete brownian_sys;
    // printf("before delete equation_systems\n");
    // delete equation_systems;
    // printf("after delete equation_systems\n");

  //delete equation_systems;
  mesh = NULL;
  pm_periodic_boundary = NULL;
  force_field = NULL;
  brownian_sys = NULL;
  //equation_systems = NULL;
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
  libmesh_example_requires(3 == LIBMESH_DIM, "--3D support");

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
EquationSystems Copss::init_system(std::string _control_fileName){
   control_fileName = _control_fileName;
   this -> read_input();
   this -> create_object_mesh();
   return this -> create_equation_systems();
}

//====================================================================
void Copss::read_input()
{
  const GetPot tmp(control_fileName);
  input_file = tmp;
  this -> read_system_info();
  this -> read_physical_info();
  this -> read_particle_info();
  this -> read_domain_info();
  this -> read_force_info();
  this -> read_ggem_info();
  this -> read_stokes_solver_info();
  this -> read_chebyshev_info();
  this -> read_run_info();
} // end read_data function

//====================================================================
void Copss::read_system_info()
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
  particle_type = input_file("particle_type", "other");
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

void Copss::read_domain_info()
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
    printf("##########################################################\n"
         "#                   Domain mesh information                      \n"
         "##########################################################\n\n");
  } // end if (comm_in.rank = 0)

  generate_mesh = input_file("generate_mesh", false); 
  if (generate_mesh){
    n_mesh.resize(input_file.vector_variable_size("n_mesh"));
    for (unsigned int i=0; i < n_mesh.size(); i++){ n_mesh[i] = input_file("n_mesh", 1, i); }
    if(comm_in.rank() == 0){
      printf(" Generate Mesh:  n_mesh = ");
      for (int i = 0 ; i < dim; i++) printf("%d, ",n_mesh[i]);
      printf("\n");
    } // end if comm_in.rank() == 0
  }
  else{
    domain_mesh_file = input_file("domain_mesh_file" , "nothing");
    if(comm_in.rank() == 0){
      printf(" Load mesh file from = %s\n", domain_mesh_file.c_str());
    } // end if comm_in.rank() == 0  
  } // end else

} // end read_domain_info()


  /*
   * read force types
   */
void Copss::read_force_info(){

  // read particle-particle force types
  num_pp_force = input_file.vector_variable_size("particle_particle_force_types");
  pp_force_type.resize(num_pp_force);
  pp_force.resize(num_pp_force);
  for (unsigned int i=0; i < num_pp_force; i++){    
    pp_force_type[i] = input_file("particle_particle_force_types", "nothing", i);
    std::vector<Real> params(input_file.vector_variable_size(pp_force_type[i]));
    if(pp_force_type[i] != "nothing"){
      for (unsigned int j = 0; j < params.size(); j++){
        params[j] = input_file(pp_force_type[i],0.0,j);
      }
    }
    pp_force[i].first = pp_force_type[i];
    pp_force[i].second = params;
  }

    // read particle-wall force types
  num_pw_force = input_file.vector_variable_size("particle_wall_force_types");
  pw_force_type.resize(num_pw_force);
  pw_force.resize(num_pw_force);
  for (unsigned int i=0; i < num_pw_force; i++){    
    pw_force_type[i] = input_file("particle_wall_force_types", "nothing" , i);
    std::vector<Real> params(input_file.vector_variable_size(pw_force_type[i]));
    if(pw_force_type[i] != "nothing"){
      for (unsigned int j = 0; j < params.size(); j++){
        params[j] = input_file(pw_force_type[i],0.0,j);
      }
    }
    pw_force[i].first = pw_force_type[i];
    pw_force[i].second = params;
  } 
  if(comm_in.rank() == 0){
    printf("##########################################################\n"
          "#    Force information (particle-particle)                \n"
          "##########################################################\n\n");
    for (int i = 0; i < num_pp_force; i++){
      printf ("  ");
      printf("%s  ", pp_force[i].first.c_str());
      for (int j = 0; j < pp_force[i].second.size(); j++){
        printf("%.6e  ", pp_force[i].second[j]);      
      }
      printf("\n");
    }
    printf("##########################################################\n"
          "#    Force information (particle-wall)                \n"
          "##########################################################\n\n");
    for (int i = 0; i < num_pw_force; i++){
      printf ("  ");
      printf("%s  ", pw_force[i].first.c_str());
      for (int j = 0; j < pw_force[i].second.size(); j++){
        printf("%.6e  ", pw_force[i].second[j]);      
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
void Copss::read_run_info(){

  //############## Without Brownian ###############################
  // For polymer_chain and bead: maximum displacement (non dimensional) of one step = max_dr_coeff * fluid mesh size minimum (hmin)
  //############## With Brownian ##################################
  // For polymer_chain: maximum displacement (non dimensional) of one step = max_dr_coeff * Ss2/Rb/Rb * 1.0
  // For bead: maximum displacement (non_dimensional) of one step = max_dr_coeff * 1.0
  with_brownian  = input_file("with_brownian", true);
  if(with_brownian){
    dt0 = input_file("dt0", 1.e-3);
    random_seed   = input_file("random_seed",111);
  }
  adaptive_dt    = input_file("adaptive_dt", true);  
  max_dr_coeff = input_file("maximum_displacement_coeff",0.1);
  restart       = input_file("restart", false);  
  restart_step  = input_file("restart_step", 0);
  restart_time  = input_file("restart_time", 0.0);
  if(restart) // update the seed for restart mode
  {
    random_seed++;
  }
  else        // set the restart_step as zero
  {
    restart_step = 0;
    restart_time = 0.0;
  }

  nstep = input_file("nstep", 1);
  write_interval = input_file("write_interval", 1);

  write_es      = input_file("write_es", true);
  out_msd_flag      = input_file("out_msd_flag", true);
  out_stretch_flag  = input_file("out_stretch_flag", false);
  out_gyration_flag = input_file("out_gyration_flag", false);
  out_com_flag      = input_file("out_com_flag", false);

  if(comm_in.rank() == 0){

    printf ("##########################################################\n"
            "#                 Run information                         \n"
            "##########################################################\n\n");
    printf("  with_brownian:  %s\n", with_brownian ? "True" : "False");
    if(with_brownian){
      printf("  random_seed = %d\n"
             "  dt0 = %.4e\n",
             random_seed, dt0);
    }
    printf("  adaptive_dt: %s\n", adaptive_dt ? "True" : "False");
    printf("  max_dr_coeff = %.4e\n", max_dr_coeff);
    printf("  write interval: %d\n", write_interval);
    printf("  Restart mode: %s\n", restart ? "True" : "False");
    printf("  restart step =  %d\n  restart time = %.4e\n", restart_step, restart_time);
    printf("  nstep = %d\n  write_interval = %d\n", nstep, write_interval);
  } // end if (comm_in.rank() == 0)

} // end read_run_info()

//============================================================================
void Copss::create_domain_mesh()
{
  if(dim == 2) {
    error_msg = "Copss::create_mesh() only works for 3D systems; 2D simulation needs extra implementation";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  mesh = new SerialMesh(comm_in);
  //mesh = std::unique_ptr<SerialMesh> (new SerialMesh (comm_in));
  if(generate_mesh and wall_type == "slit"){      
        const Real meshsize_x   = (wall_params[1] - wall_params[0])/Real( n_mesh[0] );
        const Real meshsize_y   = (wall_params[3] - wall_params[2])/Real( n_mesh[1] );
        const Real meshsize_z   = (wall_params[5] - wall_params[4])/Real( n_mesh[2] );
//        cout << "mesh_size = " <<meshsize_x  << "; "<<meshsize_y << "; "<<meshsize_z << endl;      
        min_mesh_size           = std::min(meshsize_x, meshsize_y);
        min_mesh_size           = std::min(min_mesh_size, meshsize_z);
        max_mesh_size           = std::max(meshsize_x, meshsize_y);
        max_mesh_size           = std::max(max_mesh_size, meshsize_z);
        if(comm_in.rank() == 0){
        printf("##########################################################\n"
               "#             The created mesh information                \n"
               "##########################################################\n\n"
               "   nx_mesh = %d, Lx = %.4e, hx = %.4e\n"
               "   ny_mesh = %d, Ly = %.4e, hy = %.4e\n"
               "   nz_mesh = %d, Lz = %.4e, hz = %.4e\n"
               "   minimum mesh size of fluid: hmin = %.4e\n"
               "   maximum mesh size of fliud: hmax = %.4e\n",
               n_mesh[0], wall_params[1]-wall_params[0], meshsize_x,
               n_mesh[1], wall_params[3]-wall_params[2], meshsize_y,
               n_mesh[2], wall_params[5]-wall_params[4], meshsize_z,
               min_mesh_size, max_mesh_size);
        }
     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * Create a mesh, distributed across the default MPI communicator.
     * We build a mesh with Quad9(8) elements for 2D and HEX27(20) element for 3D
    / - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        if(dim==2)
        {
        MeshTools::Generation::build_square (*mesh, n_mesh[0], n_mesh[1],
                                           wall_params[0], wall_params[1], wall_params[2], wall_params[3], QUAD8);        // QUAD8/9
        }else{
          MeshTools::Generation::build_cube (*mesh, n_mesh[0], n_mesh[1], n_mesh[2],
                                            wall_params[0], wall_params[1], wall_params[2], wall_params[3], wall_params[4], wall_params[5], HEX20);  // HEX20/27
        }
    }// end if (generate_mesh)
    else if(domain_mesh_file != "nothing"){
          mesh->read(domain_mesh_file);    
          mesh->all_second_order();
          mesh->prepare_for_use();
          const std::vector<Real> mesh_size = PMToolBox::mesh_size(*mesh);
          min_mesh_size = mesh_size[0];
          max_mesh_size = mesh_size[1];
          if(comm_in.rank() == 0){
            printf("--------------> Read finite element mesh:...\n"  
                   "##########################################################\n"
                   "#                 The Read-in mesh information            \n"
                   "##########################################################\n\n"
                   "   minimum mesh size of fluid: hmin = %.4e\n"
                   "   maximum mesh size of fliud: hmax = %.4e\n",
                   min_mesh_size, max_mesh_size);
          }
    }
    else{
          error_msg = "domain_mesh_file has to be specified";
          PMToolBox::output_message(error_msg, comm_in);
          libmesh_error();
    } 
    mesh -> print_info();
} // end function

//============================================================================
void Copss::create_periodic_boundary(){
  if(wall_type != "slit"){
    error_msg = "Copss::create_periodic_boundary() only works for wall_type = 'slit'";
    PMToolBox::output_message(error_msg, comm_in);
    libmesh_error();
  }
  const Point bbox_pmin(wall_params[0], wall_params[2], wall_params[4]);
  const Point bbox_pmax(wall_params[1], wall_params[3], wall_params[5]);
  // construct PMPeriodicBoundary class using info above
  
  pm_periodic_boundary = new PMPeriodicBoundary(bbox_pmin, bbox_pmax, periodicity, inlet, inlet_pressure);
  //pm_periodic_boundary = std::unique_ptr<PMPeriodicBoundary> 
    //                    (new PMPeriodicBoundary(bbox_pmin, bbox_pmax, periodicity, inlet, inlet_pressure));
} // end function


//=============================================================================
EquationSystems Copss::create_equation_systems()
{
  // Initialize equation_systems object using the 'mesh' we created before
  //equation_systems = new EquationSystems(*mesh);
  EquationSystems equation_systems(*mesh);
  // Add 'Stokes' system (of PMLinearImplicitSystem) to the 'equation_systems'
  PMLinearImplicitSystem& system = equation_systems.add_system<PMLinearImplicitSystem> ("Stokes");

  //Add variables to 'Stokes' system"
  u_var = system.add_variable ("u", SECOND);
  v_var = system.add_variable ("v", SECOND);
  if(dim==3)  w_var  = system.add_variable ("w", SECOND);
  const unsigned int p_var = system.add_variable ("p", FIRST);

  // attach object_mesh to pm_linear_implicit_system
  this->attach_object_mesh(system);

  // attach period boudary to pm_linear_implicit_system
  this->attach_period_boundary(system);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize the Preconditioning matrix for saddle point problems if required.
   Initialize the equation system and zero the preconditioning matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( user_defined_pc ){
    system.add_matrix("Preconditioner");
  } 
  /* Initialize the data structures for the equation system. */
  equation_systems.init();
  
  // zero the PC matrix, which MUST be done after es.init()
  if( user_defined_pc ) {
    system.get_matrix("Preconditioner").zero();
  }

  // set parameters for equation systems
  this -> set_parameters(equation_systems);

  // initialized force field
  force_field = new ForceField(system);
  system.attach_force_field(force_field);

  // print out equation system information
  output_msg = std::string("--------------> Print equation systems info\n")+
               "  System has: "+std::to_string(mesh->n_elem())+" elements,\n"+
               "              "+std::to_string(mesh->n_nodes())+" nodes,\n"+
               "              "+std::to_string(equation_systems.n_dofs())+" degrees of freedom.\n"+
               "              "+std::to_string(equation_systems.n_active_dofs())+" active degrees of freedom.\n";
  PMToolBox::output_message(output_msg, comm_in);
  equation_systems.print_info();
  return equation_systems;
}

//===============================================================================
void Copss::attach_period_boundary(PMLinearImplicitSystem& system)
{
  DofMap& dof_map = system.get_dof_map();
  /*** set PBC in x-direction ***/
  if (periodicity[0])
  {
    PeriodicBoundary pbcx(RealVectorValue(wall_params[1]-wall_params[0], 0., 0.));
    pbcx.set_variable(u_var);
    pbcx.set_variable(v_var);
    if(dim==3) pbcx.set_variable(w_var);
    //pbcx.set_variable(p_var); //*** NOT include p!
    
    // is this boundary number still true for 3D?
    if(dim==2)
    {
      pbcx.myboundary = 3;
      pbcx.pairedboundary = 1;
    }
    else if(dim==3)
    {
      pbcx.myboundary = 4; // left face (viewing from X negative to X positive)
      pbcx.pairedboundary = 2; // right face (viewing from X positive to X negative)

    } // end if
    
    dof_map.add_periodic_boundary(pbcx);
    
    PMToolBox::output_message("--->Set PBC in x direction (for 'u','v','w') finished!\n", comm_in);
    // check
    if (search_radius_p>=(wall_params[1]-wall_params[0])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in x direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Lx/2 =" + std::to_string((wall_params[1]-wall_params[0])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
    }
  }
  /*** set PBC in y-direction ***/
  if (periodicity[1])
  {
    PeriodicBoundary pbcy(RealVectorValue(0., wall_params[3]-wall_params[2], 0.));
    pbcy.set_variable(u_var);
    pbcy.set_variable(v_var);
    if(dim==3) pbcy.set_variable(w_var);
    //pbcy.set_variable(p_var); //*** NOT include p!
    
    if(dim==2)
    {
      pbcy.myboundary = 0;
      pbcy.pairedboundary = 2;
    }
    else if(dim==3)
    {
      pbcy.myboundary = 1; // bottom face, viewing from Y negative to Y positive
      pbcy.pairedboundary = 3; // top face, viewing from Y positive to Y negative
    } // end if
    
    dof_map.add_periodic_boundary(pbcy);
    // check
    if (search_radius_p>=(wall_params[3]-wall_params[2])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in y direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Ly/2 =" + std::to_string((wall_params[3]-wall_params[2])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
    }
  }
  /*** set PBC in z-direction ***/
  if (periodicity[2])
  {
    PeriodicBoundary pbcz(RealVectorValue(0., 0., wall_params[5]-wall_params[4]));
    pbcz.set_variable(u_var);
    pbcz.set_variable(v_var);
    if(dim==3) pbcz.set_variable(w_var);
    //pbcz.set_variable(p_var); //*** NOT include p!
    
    if(dim==3)
    {
      pbcz.myboundary = 0; // bottom face, viewing from Z negative to Z positive
      pbcz.pairedboundary = 5; // top face, viewing from Z positive to Z negative 
    } // end if
    
    dof_map.add_periodic_boundary(pbcz);
    // check
    if (search_radius_p>=(wall_params[5]-wall_params[4])/2.)
    {
      output_msg = std::string("****************************** warning: ********************************\n")+
                   "**** The search radius is larger than half domain length in z direction\n!"+
                   "**** search radius = "+ std::to_string(search_radius_p) + 
                   ", half domain size Lx/2 =" + std::to_string((wall_params[5]-wall_params[4])/2.)+"\n"+
                   "************************************************************************\n\n";
      PMToolBox::output_message(output_msg, comm_in);
    }
  }  
}


//============================================================================================
void Copss::solve_undisturbed_system(EquationSystems& equation_systems)
{
   // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute undisturbed velocity field without particles.
   NOTE: We MUST re-init particle-mesh before solving Stokes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  system.reinit_system();
  reinit_stokes = true;
  system.solve_stokes("undisturbed",reinit_stokes);
  v0_ptr = system.solution->clone(); // backup v0
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write out the equation systems if write_es = true at Step 0 (undisturbed field)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  exodus_ptr = new ExodusII_IO(*mesh);
  if(write_es && restart==false)
  {
    //system.add_local_solution(); // Don't add local solution for undisturbed system!
#ifdef LIBMESH_HAVE_EXODUS_API
    exodus_ptr->write_equation_systems(out_system_filename+".e", equation_systems);
#endif
  }
}

//============================================================================================
void Copss::create_brownian_system(EquationSystems& equation_systems)
{
  brownian_sys = new BrownianSystem(equation_systems);
  brownian_sys->init_petsc_random(&rand_ctx);
  brownian_sys->_create_shell_mat(n_vec, &M);
  brownian_sys->_create_petsc_vec(n_vec,&R0);
  VecDuplicate(R0,&U0);
  VecDuplicate(R0,&R_mid);
  VecDuplicate(R0,&dw_mid);
  brownian_sys->extract_particle_vector(&ROUT,"coordinate","extract");
  VecDuplicate(ROUT,&RIN);
  VecCopy(ROUT,RIN);  // RIN = ROUT = the initial position vector
  brownian_sys->set_std_random_seed(random_seed);  // random seed
 
  if(restart)
  {
    // read RIN & ROUT from local file output during the previous simulation
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_RIN.dat",FILE_MODE_READ,&viewer);
    VecLoad(RIN,viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_READ,&viewer);
    VecLoad(ROUT,viewer);
  }
  else
  {
    // write out binary file of RIN, which may be used at restart mode.
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_RIN.dat",FILE_MODE_WRITE,&viewer);
    VecView(RIN,viewer);
  }
  comm_in.barrier();
  center0 = brownian_sys->center_of_mass(RIN);
  if(!restart)
  {
    /* Output mean square displacement and radius of gyration at step 0 (the origin) */
    brownian_sys->output_statistics_step0(out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag, RIN);
  }  
}

//==============================================================================================
void Copss::fixman_integrate(EquationSystems& equation_systems, unsigned int i)
{

  // get stokes system from equation systems
  PMLinearImplicitSystem& system = equation_systems.get_system<PMLinearImplicitSystem> ("Stokes");
 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the "disturbed" particle velocity + "undisturbed" velocity = U0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //cout <<"Compute the disturbed particle velocity at step "<<i+1<<endl;
    if(i>0){
      *(system.solution) = *v0_ptr; // re-assign the undisturbed solution
      // Update the local values to reflect the solution on neighboring processors
      system.update();
      system.reinit_system();
    }
    // compute undisturbed velocity of points 
    system.compute_point_velocity("undisturbed", vel0);
    reinit_stokes = false;
    system.solve_stokes("disturbed",reinit_stokes); // Using StokesSolver
    // compute distrubed velocity of points
    system.compute_point_velocity("disturbed", vel1);
    // add up undistrubed and disturbed velocity of points
    for(std::size_t j=0; j<vel1.size();++j) vel1[j] += vel0[j];
    // transform total point velocity to U0 in Brownian_system
    brownian_sys->vector_transform(vel1, &U0, "forward");
 
    /*---------------------------------------------------------------------------------------
     * write equation system at step i
    -----------------------------------------------------------------------------------------*/
    if(i%write_interval == 0){
      o_step++;
      if( write_es){
        system.add_local_solution();  // add local solution for the disturbed system
        system.solution->add(*v0_ptr);// add the undisturbed solution
#ifdef LIBMESH_HAVE_EXODUS_API
        exodus_ptr->append(true);
        exodus_ptr->write_timestep(out_system_filename+".e",system.get_equation_systems(),o_step,o_step);
#endif
      } // end if (write es)   
    } // end if (i % write_interval == 0 )
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adaptive time step.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Real dt = 0;
    hmin = equation_systems.parameters.get<Real>("fluid mesh size");
    if(adaptive_dt){
      Real vp_max = 0.0, vp_min = 0.0;
      for(unsigned int k=0; k<dim;++k) {
        vp_max += vel1[k]*vel1[k];
      }
      vp_min = vp_max;
      for(unsigned int j=1; j<NP;++j)
      {
        Real vp_norm = 0.0;
        for(std::size_t k=0; k<dim;++k) vp_norm += vel1[j*dim+k]*vel1[j*dim+k];
        vp_max = std::max(vp_max,vp_norm);
        vp_min = std::min(vp_min,vp_norm);
      }
      vp_max = std::sqrt(vp_max);     // maximum magnitude of particle velocity
      vp_min = std::sqrt(vp_min);
      if(with_brownian) {
        dt = (vp_max == 0) ? (max_dr_coeff) : (max_dr_coeff * 1. / vp_max);  // maximum |dr| =  max_dr = dt0 * 1; dt0 = 0.1 for beads; dt0 = 0.1*Ss2/Rb/Rb for polymers
      }
      else {
        dt = (vp_max == 0) ? (max_dr_coeff * hmin) : (max_dr_coeff * hmin / vp_max); // maximum |dr| = dt0 * hmin
      }
      if(i % write_interval == 0){
        output_msg=
                "# Max velocity magnitude : "+ std::to_string(vp_max) +"\n"+
                "# Min velocity magnitude : "+ std::to_string(vp_min) +"\n"+
                "# minimum fluid mesh size :"+ std::to_string(hmin) +"\n"+
                "# The adaptive time increment at step "+std::to_string(i+1) + " dt : "+ std::to_string(dt);
        PMToolBox::output_message(output_msg, comm_in);
      } // end if (i% write_interval)  
    }
    else{
      if(with_brownian) {
        dt = max_dr_coeff * 1;  // maximum |dr| =  max_dr = dt0 * 1; dt0 = 0.1 for beads; dt0 = 0.1*Ss2/Rb/Rb for polymers
      }
      else {
        dt = max_dr_coeff * hmin; // maximum |dr| = dt0 * hmin
      } // end if (with brownian)

      if(i % write_interval == 0){
        output_msg = "# The fixed time increment at step "+std::to_string(i+1) +" dt : " +std::to_string(dt)+"\n";
        PMToolBox::output_message(output_msg,comm_in);
      } // end if (i % write_interval == 0)
    } // end if (adaptive_dt)

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      If with Brownian motion, we use midpoint scheme
      If without Brownian motion, we use normal stepping: dR = Utotal*dt
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */  
    if (with_brownian)
    {
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Generate random vector dw whose mean = 0, variance = sqrt(2*dt)
       petsc_random_vector generates a uniform distribution [0 1] whose
       mean = 0.5 and variance = 1/12, so we need a shift and scale operation.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      Real mean_dw = 0.0, variance_dw = 0.0;
      
      // A more precise way is to construct a random vector with gaussian distribution
      const Real std_dev  = std::sqrt(dt);
      brownian_sys->std_random_vector(0.0,std_dev,"gaussian",&dw);
      brownian_sys->_vector_mean_variance(dw, mean_dw, variance_dw);
      VecScale(dw,std::sqrt(2.0));
      
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Print out the mean and variance or view the generated vector.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      //PetscPrintf(PETSC_COMM_WORLD,
      //            "Generated random_vector:               mean = %f, variance = %f\n",
      //            mean_dw, variance_dw);
      //PetscPrintf(PETSC_COMM_WORLD,
      //            "Exact values for uniform distribution: mean = %f, variance = %f\n",
      //            0.5, 1./12.);
      
      
      // Compute dw = B^-1 * dw using Chebyshev polynomial, dw will be changed!
      VecCopy (dw,dw_mid);  // save dw to dw_mid, which will be used for Chebyshev
      for(std::size_t j=0; j<2; j++)
      {
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the max/min eigenvalues if needed. Otherwise, magnify the interval.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        if(compute_eigen){
          //cout << "Compute the max & min eigenvalues for Chebyshev polynomial at step "<<i+1<<endl;
          brownian_sys->compute_eigenvalues(eig_min,eig_max,tol_eigen);
        }
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the Brownian displacement B^-1 * dw using Chebyshev approximation.
         Here dw is both input and output variables, so it will be changed.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        //PetscPrintf(PETSC_COMM_WORLD,
        //            "--->eig_min = %f, eig_max = %f, tol_cheb = %f, max_n_cheb = %d\n",
        //            eig_min,eig_max,tol_cheb,max_n_cheb);
        cheb_converge = brownian_sys->chebyshev_polynomial_approximation(max_n_cheb,
                                                                        eig_min,eig_max,tol_cheb,&dw);       
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         If converged, dw returns the Brownian displacement, then break the j-loop;
         Otherwise, recompute eigenvalues
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        if(cheb_converge){
          compute_eigen = false; break;
        }
        else{
          compute_eigen = true;
          VecCopy(dw_mid,dw); /*copy back, recompute eigenvalues*/
          //cout << "It is necessry to re-compute the eigenvalues at step " <<i+1<<endl;
        }
      } // end for j-loop

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Double-check the convergence of Chebyshev polynomial approximation
       *** If cheb does NOT converge, consider re-generating a rand vector!
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      if(!cheb_converge)
      {
        PMToolBox::output_message("****** Warning: Chebysheve failed to converge at step " +std::to_string(i+1), comm_in);
      }

      // Magnify the spectral range by a factor (1.05 by default).
      eig_max *= eig_factor; eig_min /= eig_factor;
      
      // Compute dw_mid = D*B^-1*dw, which can be obtained by solving the Stokes
      brownian_sys->hi_ewald(M,dw,dw_mid);  // dw_mid = D * dw

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Particle coordinate vector R0.
       Move the particle R_mid = R0 + 0.5*(U0+U1)*dt (deterministic)
       and R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw     (stochastic)
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
      VecWAXPY(R_mid,0.5*dt,U0,R0);  // R_mid = R0 + 0.5*dt*(U0+U1)  (R0 and U0 do NOT change)
      coef = 0.5;                    // coefficient. sqrt(2) is introduced when generating dw
      VecAXPY(R_mid,coef,dw_mid);    // R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw
      brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords
   
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Check and correct beads' position at the midpoint
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      force_field->check_walls(); // check pbc and inpenetrable wall
      this -> update_object("after midpoint at step"+std::to_string(i+1));
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Update the particle mesh for the mid-point step,
       and recompute U0 + U1_mid, D_mid*(B^-1*dw)
       NOTE: the FEM solution of undisturbed field doesn't change, but particles
       move, so U0 needs to be re-evaluated at the new position.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      *(system.solution) = *v0_ptr;       // re-assign the undisturbed solution
      system.update();
      system.reinit_system();
      system.compute_point_velocity("undisturbed", vel0);
      reinit_stokes = false;
      system.solve_stokes("disturbed",reinit_stokes);   // solve the disturbed solution
      system.compute_point_velocity("disturbed", vel1);
      for(std::size_t j=0; j<vel1.size();++j) vel1[j] += vel0[j];
      brownian_sys->vector_transform(vel1, &U0, "forward"); // (U0+U1)_mid
      brownian_sys->hi_ewald(M,dw,dw_mid);  // dw_mid = D_mid*dw, where dw=B^-1*dw computed above
   
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       the mid-point to the NEW point, and update the particle coordinates
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      //cout << "Update from mid-point to the NEW particle coordinates at step " <<i+1<<endl;
      VecWAXPY(R_mid,dt,U0,R0);         // R_mid = R0 + dt*U0_mid
      VecAXPY(R_mid,2.0*coef,dw_mid); // R_mid = R_mid + sqrt(2)*D_mid*B^-1*dw
      brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign");
  
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Check and correct the beads' position again after the midpoint update
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      force_field->check_walls(); // check pbc and inpenetrable wall
      this -> update_object("after step "+std::to_string(i+1));

      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT,dt,U0);            // ROUT = ROUT + dt*U0_mid
      VecAXPY(ROUT,2.0*coef,dw_mid);  // ROUT = ROUT + sqrt(2)*D_mid*B^-1*dw
    } // end if Brownian
    else{ // if without Brownian

      // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
      brownian_sys->extract_particle_vector(&R0,"coordinate","extract");
      VecWAXPY(R_mid,dt,U0,R0);  // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
      brownian_sys->extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords
 
      // Check and correct beads' position at the midpoint
      force_field->check_walls(); // check pbc and inpenetrable wall
      this -> update_object("after step "+std::to_string(i+1));

      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT,dt,U0); // ROUT = ROUT + dt*U0_mid   
    } // end else (without_brownian)
    real_time += dt;
}


void Copss::destroy()
{
  MatDestroy(&M);
  VecDestroy(&U0);
  VecDestroy(&R0);
  VecDestroy(&R_mid);
  VecDestroy(&dw_mid);
  PetscRandomDestroy(&rand_ctx);
  if(with_brownian){
    VecDestroy(&dw);
  }
  if(exodus_ptr) {
    delete exodus_ptr;
  }
  PetscViewerDestroy(&viewer);
}

} // end namespace





