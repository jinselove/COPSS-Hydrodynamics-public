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



// C++ includes
#include <iostream>
#include <algorithm>
#include <cstring>
#include <math.h>
#include <time.h>
#include <fstream> 
#include <tuple>
#include <stdlib.h>

// Libmesh includes 
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/getpot.h"
#include "libmesh/exodusII_io.h"

#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/dof_map.h"
#include "libmesh/linear_solver.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"

#include "libmesh/slepc_macro.h"

// User defined includes
#include "../../src/point_particle.h"
#include "../../src/particle_mesh.h"
#include "../../src/point_mesh.h"
#include "../../src/force_field.h"
#include "../../src/pm_linear_implicit_system.h"
#include "../../src/brownian_system.h"
#include "../../src/pm_periodic_boundary.h"
#include "../../src/chebyshev.h"
#include "../../src/pm_toolbox.h"
#include "../../src/polymer_chain.h"
#include "../../src/random_generator.h"
#include "../../src/stokes_solver.h"
#include "../../src/ggem_system.h"


/*! This file serves as a template to build\solve a Stokes system
    using our pFE-GgEm code 
*/

// Bring in everything from the libMesh namespace
using namespace libMesh;
using std::cout;
using std::endl;
using std::string;

// ----------------------------- check requires for libmesh  -----------------------------
int check_libmesh()
{
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

/*******************************************************************************************
 *				The main program					   *
 *******************************************************************************************/
int main (int argc, char** argv)
{
  
  //============================0.  Initialize libMesh.===========================
  LibMeshInit init (argc, argv);
  
  // Check libMesh requirements 
  check_libmesh();
  
  // Print out the current date and time 
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  const Parallel::Communicator& comm_in = init.comm();
  cout<<"---------------------------------------------------------------------"<<endl;
  cout<<"The program starts."<<endl;
  cout<<"The Current data/time is: "<<asctime(timeinfo);
  cout<<"------------------------------------------------------------------------"<<endl<<endl;
  
  cout<<"============================0.  Initialize libMesh.==========================="<<endl<<endl;
  

//============================1. Read input parameters ============================
  cout<<"============================1. read input parameters ============================"<<endl<<endl;

//===========1/1 Read polymer_control file and define simulation parameters
  cout<<"==>(1/1) Read polymer_control file and define simulation parameters"<<endl<<endl;
  GetPot input_file("polymer_control.in");
  const std::string test_name = input_file("test_name", "validation");	
  const Real kB = 1.380662E-17;//1.380662E-23(J/K) = 1.3806623E-23(N*m/K) = 1.380662E-17(N*um/K)
  const Real T = input_file("temperature", 297);// K
  const Real PI     = libMesh::pi;
  const Real kBT    = kB * T; //(N*um)  

  // (1) Stokes solver parameters
  const int max_linear_iterations = input_file("max_linear_iterations", 100);
  const Real linear_solver_rtol   = input_file("linear_solver_rtol", 1E-6);
  const Real linear_solver_atol   = input_file("linear_solver_atol", 1E-6);
  bool user_defined_pc            = input_file("user_defined_pc", true);
  const bool schur_user_ksp       = input_file("schur_user_ksp", false);
  const Real schur_user_ksp_rtol  = input_file("schur_user_ksp_rtol", 1E-6);
  const Real schur_user_ksp_atol  = input_file("schur_user_ksp_atol", 1E-6);
  const std::string schur_pc_type = input_file("schur_pc_type", "SMp");
  const std::string stokes_solver_type = input_file("stokes_solver", "superLU_dist");
  StokesSolverType solver_type;
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
  
  // (2) Geometry parameters
  const unsigned int dim          = input_file("dimension", 3);
  //=============== wall type and radius of spherical cavity
  const string wall_type = input_file("wall_type", "not_defined");
  std::vector<Real> wall_params(input_file.vector_variable_size(wall_type));
  if(wall_type != "not_defined"){
    for (unsigned int j = 0; j < wall_params.size(); j++){
      wall_params[j] = input_file(wall_type,0.0,j);
    }
  }
  else{
    cout <<endl<<"-----------------------------------------------------------------------"<<endl
               <<"warning: check the wall_type definition in control file (1. wall_type; 2. wall_params)"<<endl
               <<"------------------------------------------------------------------------"<<endl;
    libmesh_error();
  }
  std::vector<Real> XYZAB(dim*2);
  if(wall_type == "slit") XYZAB = wall_params;

  //=============== alpha
  Real alpha                = input_file("alpha", 0.1);
  //=============== periodicity
  std::vector<bool> periodicity;
  periodicity.resize(input_file.vector_variable_size("periodicity"));
  for (unsigned int i=0; i < periodicity.size(); i++){ periodicity[i] = input_file("periodicity", false, i); }
  if(periodicity[0]==true and periodicity[1] == true and periodicity[2]==true){
  	cout <<endl<<"-----------------------------------------------------------------------"<<endl
  	           <<"warning: The box cannot be periodic on all directions at the same time. (required by FEM)"<<endl
  	           <<"------------------------------------------------------------------------"<<endl;
  	libmesh_error();
  }
  //============== inlet 
  std::vector<bool> inlet;
  inlet.resize(input_file.vector_variable_size("inlet"));
  for (unsigned int i=0; i < inlet.size(); i++){ 
     inlet[i] = input_file("inlet", false, i);
     if(inlet[i]==true and periodicity[i]==true) {
	    cout <<endl<<"-----------------------------------------------------------------------"<<endl
	           <<"warning: A inlet direction has to be non-periodicity"<<endl
	           <<"------------------------------------------------------------------------"<<endl;
	    libmesh_error();
     }
  }
  //============== inlet pressure
  std::vector<Real> inlet_pressure;
  inlet_pressure.resize(input_file.vector_variable_size("inlet_pressure"));
  for (unsigned int i=0; i < inlet_pressure.size(); i++){ inlet_pressure[i] = input_file("inlet_pressure", 0, i); }

  // (3) Physical parameters of polymer and fluids
  
  // Universal parameters for all particle models (point_particle, polymer_chain, finite_size_particle, etc.)
  const Real viscosity            = input_file("viscosity", 1.0); // viscosity (cP = N*s/um^2)
  const Real Rb                   = input_file("radius", 0.10); // radius of the bead (um)
  const Real  drag_c      = 6.*PI*viscosity*Rb;    // Drag coefficient (N*s/um)
  const Real  Db          = kBT/drag_c;     // diffusivity of a bead (um^2/s)
 
  // Specific parameters
  unsigned int Nb; //for "polymer_chain" and "point_particle"
  unsigned int Ns, nChains; // for "polymer_chain"
  Real bk,Nks,q0,chain_length,Dc, Ss2;// for "polymer_chain" models
 
  // read parameters depending on particle_model (point_particle only)
  const string particle_type = input_file("particle_type","other");
  if(particle_type != "point_particle"){
    cout << endl << "==================================================================" << endl
         <<"Invalid particle type (" <<particle_type <<" )"<<endl
         <<"=================================================================="<<endl;
    libmesh_error();
  }

  const string point_particle_model	= input_file("point_particle_model", "other");
  if(point_particle_model == "bead"){
  	Nb = input_file("Nb",21); // total # of point particles
  	Ns = Nb - 1;//we can look at these beads as a single chain without spring force. Just for development convinence
  }
  else if(point_particle_model == "polymer_chain"){
    Nb = input_file("Nb", 21);// total # of beads
  	Ns = input_file("Ns", 20);// # of springs per Chain
  	nChains = Nb / (Ns+1);
  	bk = input_file("bk", 1E-6);//Kuhn length(um)
  	Nks = input_file("Nks", 1E-6);// # of Kuhn length per spring
  	q0 = Nks * bk;// maximum spring length (um)
  	chain_length = Ns * q0; // contour length of the spring (um)
  	Dc = Db / Real(Nb); // Diffusivity of the chain (um^2/s)
  	Ss2 = Nks*bk*bk/6.; // (um^2)
  }
  else{
	  cout << endl << "==================================================================" << endl
	     <<"Invalid point particle type"<<endl
	     <<"=================================================================="<<endl;
    libmesh_error();
  }
  //make sure we have the right combination of Nb and Ns
  if((Nb % (Ns + 1)) != 0){
    cout << endl << "==================================================================" << endl
       <<"Incorrect combination of Nb (number of beads) and Ns (number of springs per chain)"<<endl
       <<"=================================================================="<<endl;
    libmesh_error();     
  }  

// read force field to the system
  // read particle-particle force types
  unsigned int num_pp_force = input_file.vector_variable_size("particle_particle_force_types");
  std::vector<string> pp_force_types(num_pp_force);
  // use the type "ForceField::type_force" defined in force_field.h
  std::vector<ForceField::type_force> pp_force(num_pp_force);
  for (unsigned int i=0; i < num_pp_force; i++){    
      pp_force_types[i] = input_file("particle_particle_force_types", "nothing", i);
      std::vector<Real> params(input_file.vector_variable_size(pp_force_types[i]));
      if(pp_force_types[i] != "nothing"){
        for (unsigned int j = 0; j < params.size(); j++){
          params[j] = input_file(pp_force_types[i],0.0,j);
        }
      }
      pp_force[i].first = pp_force_types[i];
      pp_force[i].second = params;
  }

  // read particle-wall force types
  unsigned int num_pw_force = input_file.vector_variable_size("particle_wall_force_types");
  std::vector<string> pw_force_types(num_pw_force);
  std::vector<ForceField::type_force> pw_force(num_pw_force);
  for (unsigned int i=0; i < num_pw_force; i++){    
      pw_force_types[i] = input_file("particle_wall_force_types", "nothing" , i);
      std::vector<Real> params(input_file.vector_variable_size(pw_force_types[i]));
      if(pw_force_types[i] != "nothing"){
        for (unsigned int j = 0; j < params.size(); j++){
          params[j] = input_file(pw_force_types[i],0.0,j);
        }
      }
      pw_force[i].first = pw_force_types[i];
      pw_force[i].second = params;
  } 

  // (4) characteristic variables
  const Real tc   = drag_c*Rb*Rb/kBT;       // diffusion time (s)
  const Real uc   = kBT/(drag_c*Rb);        // characteristic velocity (um/s)
  const Real fc   = kBT/Rb;                 // characteristic force (N)
  const Real muc  = 1./(6.*PI);             // non-dimensional viscosity

  // (5) parameters for dynamics process
  const unsigned int nstep  = input_file("nstep", 100);

  //############## Without Brownian ###############################
  // For polymer_chain and bead: maximum displacement (non dimensional) of one step = max_dr_coeff * fluid mesh size minimum (hmin)
  //############## With Brownian ##################################
  // For polymer_chain: maximum displacement (non dimensional) of one step = max_dr_coeff * Ss2/Rb/Rb * 1.0
  // For bead: maximum displacement (non_dimensional) of one step = max_dr_coeff * 1.0
  Real max_dr_coeff = input_file("maximum_displacement_coeff",0.1);
  if(point_particle_model == "polymer_chain") max_dr_coeff *= Ss2/Rb/Rb; 
  const bool adaptive_dt    = input_file("adaptive_dt", true);
  const unsigned int write_interval = input_file("write_interval", 1);
  const bool  restart       = input_file("restart", false);
  std::size_t restart_step  = input_file("restart_step", 0);
  std::size_t random_seed   = input_file("random_seed",111);
  Real        restart_time  = input_file("restart_time", 0.0);
  if(restart) // update the seed for restart mode
  {
    random_seed++;
  }
  else        // set the restart_step as zero
  {
    restart_step = 0;
    restart_time = 0.0;
  }

  const bool with_brownian  = input_file("with_brownian", true);
  bool compute_eigen;
  unsigned int max_n_cheb;
  Real tol_cheb, tol_eigen, eig_factor;
  if(with_brownian){
    // by default, compute_eigen = true since we need to compute eigen values at step 1
    compute_eigen = true;
    max_n_cheb = input_file("max_n_cheb", 10);
    tol_cheb = input_file("tol_cheb", 0.1);
    eig_factor = input_file("eig_factor", 1.05);
    tol_eigen = input_file("tol_eigen", 0.01);
  }
  const bool  write_es      = input_file("write_es", true);
  const bool  print_info    = input_file("print_info", false);
  bool    out_msd_flag      = input_file("out_msd_flag", true);
  bool    out_stretch_flag  = input_file("out_stretch_flag", false);
  bool    out_gyration_flag = input_file("out_gyration_flag", false);
  bool    out_com_flag      = input_file("out_com_flag", false);
  std::ostringstream oss;


  cout << "##########################################################"<<endl
       << "#                       system_name                       "<<endl
       << "##########################################################"<<endl<<endl;
	cout << "-----------> system_name: "<<test_name << endl;
  
  cout <<endl<< "##########################################################"<<endl
       << "#                  Physical Parameters                    "<<endl
       << "##########################################################"<<endl<<endl;
  // for all particle models
  cout << "   temperature           T   = " << T << "(K)" <<endl;
  cout << "   viscosity             mu  = " << viscosity <<" (cP = N*s/um^2)"<<endl;
  cout << "   Particle type             : "<<particle_type <<endl;
  cout << "   Point Particle model      : "<<point_particle_model <<endl;
  cout << "   Energy unit           kBT = " << kBT << " (N*um = N*um)"<<endl;
  cout << "   Radius of the bead     a  = " << Rb << " (um)"<<endl;
  cout << "   bead diffusivity      Db  = " << Db << " (um^2/s)"<<endl;
  cout << "   HI Drag coefficient  zeta = 6*PI*mu*a = " << drag_c << " (N*s/um)"<<endl;
  cout << "   ksi = sqrt(PI)/(3a)       = " << std::sqrt(PI)/(3.*Rb) <<"(1/um)" <<endl;
  cout << "   number of beads       Nb  = " << Nb << endl;

  // for particular models
  if(point_particle_model == "polymer_chain"){
      cout << "   number of springs per Chain         Ns  = " << Ns << endl;
        cout << "   number of Chains              nChains = " <<nChains <<endl;
        cout << "   Kuhn length                       bk  = " << bk << " (um)"<<endl;
        cout << "   # of Kuhn segment per spring      Nks = " << Nks << endl;
        cout << "   second moment of polymer chain    Ss2 = " << Ss2 << " (um^2)"<<endl;
        cout << "   maximum spring length             q0  = " << q0 << " (um)"<<endl;
        cout << "   chain length of polymer           Lc  = " << chain_length << " (um)"<<endl;
        cout << "   chain diffusivity                 Dc  = " << Dc << " (um^2/s)" << endl;
  }

    cout << endl<<"------------> The characteristic variables:"<<endl;
    cout << "   characteristic time          = " << tc << " (s)"<<endl;
    cout << "   characteristic velocity      = " << uc << " (um/s)"<<endl;
    cout << "   characteristic force         = " << fc << " (N)"<<endl;


    cout << endl<<"------------> The non-dimensional variables:"<<endl;

  if(point_particle_model == "polymer_chain"){
      cout << "   non-dimensional Kuhn length    bk/a     = " << bk/Rb <<endl;
      cout << "   non-dimensional spring length  q0/a     = " << q0/Rb <<endl;
      cout << "   non-dimensional contour length Lc/a     = " << chain_length/Rb <<endl;
      cout << "   non-dimensional Ss/a = sqrt(Ss2/a^2)    = " << std::sqrt(Ss2/Rb/Rb) <<endl;
      cout << "   non-dimensional ksi = sqrt(PI)/(3a0)    = " << std::sqrt(PI)/(3.) <<endl;
  }
      cout << "   non-dimensional bead radius      a0     = " << 1.0 << endl;
      cout << "   non-dimensional ksi = sqrt(PI)/(3a0)    = " << std::sqrt(PI)/(3.) <<endl;
      cout <<endl<<endl;

  cout <<endl<< "##########################################################"<<endl
       << "#                  Geometry information                   " <<endl
       << "##########################################################"<<endl<<endl;
  cout << "-----------> Dimension: " << dim << endl;
  cout << "-----------> Wall type: " << wall_type << endl;
  cout << "-----------> Wall size parameters: ";
       for (int i = 0; i < wall_params.size(); ++i)
       {
         cout << wall_params[i] <<"   ";
       }
       cout << endl;
  cout << "-----------> Periodicity of the box: "<< std::boolalpha << periodicity[0] << ", "<<std::boolalpha << periodicity[1]<<", "<<std::boolalpha <<periodicity[2]<<endl
       << "-----------> Inlet/Outlet of the box: "<< std::boolalpha << inlet[0] <<"(pressure = " <<inlet_pressure[0] <<" ), "
      << std::boolalpha << inlet[1] <<"(pressure = " <<inlet_pressure[1] <<" ), " 
      << std::boolalpha << inlet[2] <<"(pressure = " <<inlet_pressure[2] <<" )" << endl;
  cout <<endl<< "##########################################################"<<endl
       << "#                   Mesh information                     " <<endl
       << "##########################################################"<<endl<<endl;
  bool generate_mesh = input_file("generate_mesh", true); 
  const string domain_mesh_file = input_file("domain_mesh_file","nothing");
  std::vector<Real> n_mesh;
  cout <<"-----------> generate_mesh (cubit mesh)="<<std::boolalpha<<generate_mesh << endl;
  if(generate_mesh){
    n_mesh.resize(input_file.vector_variable_size("n_mesh"));
    for (unsigned int i=0; i < n_mesh.size(); i++){ n_mesh[i] = input_file("n_mesh", 10, i); }
    cout<<"-----------> n_mesh = "<<n_mesh[0] << "; "<< n_mesh[1] << "; "<<n_mesh[2] << endl;
  }//end if generate_mesh
  else{
   cout<<"-----------> Read mesh file from "<<domain_mesh_file<<endl;
  }
  cout <<endl<< "##########################################################"<<endl
       << "#    Force information (particle-particle)                     " <<endl
       << "##########################################################"<<endl<<endl;
  for (int i = 0; i < num_pp_force; i++){
    cout << "-----------> ";
    cout << pp_force[i].first <<"= '";
    for (int j = 0; j < pp_force[i].second.size(); j++){
          cout << pp_force[i].second[j] <<"  ";       
    }
    cout << "'" <<endl;
  }

  cout <<endl<< "##########################################################"<<endl
       << "#    Force information (particle-wall)                     " <<endl
       << "##########################################################"<<endl<<endl;
  for (int i = 0; i < num_pw_force; i++){
    cout << "-----------> ";
    cout << pw_force[i].first <<"= '";
    for (int j = 0; j < pw_force[i].second.size(); j++){
          cout << pw_force[i].second[j] <<"  ";       
    }
    cout << "'" <<endl;
  }

  cout << endl<<"##########################################################"<<endl
       << "#                 GGEM information                      " <<endl
       << "##########################################################"<<endl<<endl;
  
  cout << "-----------> the smoothing parameter in GGEM alpha = " << alpha << endl; 
  cout << "-----------> recommend meshsize <= " << 1./(std::sqrt(2)*alpha) <<endl;

  cout <<endl<< "##########################################################"<<endl
       << "#                 Solver information                      " <<endl
       << "##########################################################"<<endl<<endl;
  cout << "-----------> Stokes solver type = " << stokes_solver_type <<endl;
  if (stokes_solver_type=="field_split"){
    cout<<"-----------> FieldSplit Schur Complement Reduction Solver"<<endl;
    cout<<"-----------> schur_pc_type = " << schur_pc_type << endl;
      if(schur_user_ksp){
            cout<<"----------->  user defined KSP is used for Schur Complement!"<< endl;
            cout<<"----------->  KSP rel tolerance for Schur Complement solver is = " << schur_user_ksp_rtol <<endl;
            cout<<"----------->  KSP abs tolerance for Schur Complement solver is = " << schur_user_ksp_atol <<endl;
        }// end if(schur_user_ksp)
    }// end if (stokes_solver_type == "field_split")

  if(with_brownian){
    cout <<endl<< "##########################################################"<<endl
         << "#   Chebyshev information (only for brownian System)            " <<endl
         << "##########################################################"<<endl<<endl;  
    cout << "-----------> compute eigen values  = " <<std::boolalpha << compute_eigen <<endl;
    cout << "-----------> max number of chebyshev polynomial = " <<max_n_cheb <<endl;
    cout << "-----------> tolerance of chebyshev polynomial = " <<tol_cheb <<endl;
    cout << "-----------> factor of eigenvalues range = " <<eig_factor <<endl;
    cout << "-----------> tolerance of eigenvalues convergence = " <<tol_eigen <<endl;
  }

  cout <<endl<< "##########################################################"<<endl
       << "#                 Run information                      " <<endl
       << "##########################################################"<<endl<<endl;
  cout << "-----------> with_brownian: " <<std::boolalpha<<with_brownian <<endl;
  cout << "-----------> adaptive_dt: " << std::boolalpha << adaptive_dt << endl;
  cout << "-----------> max_dr_coeff: " << max_dr_coeff << endl;
  cout << "-----------> write interval: " <<write_interval <<endl;
  cout << "-----------> Restart mode: "<<std::boolalpha << restart <<"; restart step: "<<restart_step <<"; restart time: "<<restart_time <<endl;
  cout << "-----------> random seed: " <<random_seed <<endl;
  cout << "-----------> nstep = " <<nstep <<endl;

   
   
   
  cout<<endl<<"============================2. Create Point-mesh object ============================"<<endl<<endl;
  //===============> (1/4) Generate/Create Mesh
  cout<<"==>(1/4) Generate/Create Mesh"<<endl;
  SerialMesh mesh(comm_in);   //Mesh mesh(comm_in);

  Real min_mesh_size, max_mesh_size;  
  if(generate_mesh){      
      const Real meshsize_x   = (XYZAB[1] - XYZAB[0])/Real( n_mesh[0] );
      const Real meshsize_y   = (XYZAB[3] - XYZAB[2])/Real( n_mesh[1] );
      const Real meshsize_z   = (XYZAB[5] - XYZAB[4])/Real( n_mesh[2] );
      cout << "mesh_size = " <<meshsize_x  << "; "<<meshsize_y << "; "<<meshsize_z << endl;      
      min_mesh_size           = std::min(meshsize_x, meshsize_y);
      min_mesh_size           = std::min(min_mesh_size, meshsize_z);
      max_mesh_size           = std::max(meshsize_x, meshsize_y);
      max_mesh_size           = std::max(max_mesh_size, meshsize_z);
      cout <<endl<< "##########################################################"<<endl
       << "#                 The created mesh information                      " <<endl
       << "##########################################################"<<endl<<endl;
      cout << "   nx_mesh = " << n_mesh[0] <<", Lx = " << XYZAB[1]-XYZAB[0] <<", hx = "<< meshsize_x <<endl
           << "   ny_mesh = " << n_mesh[1] <<", Ly = " << XYZAB[3]-XYZAB[2] <<", hy = "<< meshsize_y <<endl
           << "   nz_mesh = " << n_mesh[2] <<", Lz = " << XYZAB[5]-XYZAB[4] <<", hz = "<< meshsize_z <<endl
           << "   minimum mesh size of fluid: hmin = " << min_mesh_size << endl
           << "   maximum mesh size of fluid: hmax = " << max_mesh_size <<endl;
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Create a mesh, distributed across the default MPI communicator.
   * We build a mesh with Quad9(8) elements for 2D and HEX27(20) element for 3D
  / - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 	    cout << "--------------> Create finite element mesh..."<<endl;
    	if(dim==2)
    	{
  		MeshTools::Generation::build_square (mesh, n_mesh[0], n_mesh[1],
                                         XYZAB[0], XYZAB[1], XYZAB[2], XYZAB[3], QUAD8);        // QUAD8/9
    	}else if(dim==3)
    	{
  			MeshTools::Generation::build_cube (mesh, n_mesh[0], n_mesh[1], n_mesh[2],
                                       XYZAB[0], XYZAB[1], XYZAB[2], XYZAB[3], XYZAB[4], XYZAB[5], HEX20);  // HEX20/27
    	}else{libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");}
  }// end if (generate_mesh)
  else{
      if(domain_mesh_file != "nothing"){
        mesh.read(domain_mesh_file);
        cout << "--------------> Read finite element mesh:..."<<endl;      
        mesh.all_second_order();
        mesh.prepare_for_use();
        const std::vector<Real> mesh_size = PMToolBox::mesh_size(mesh);
        min_mesh_size = mesh_size[0];
        max_mesh_size = mesh_size[1];
      cout <<endl<< "##########################################################"<<endl
       << "#                 The Read-in mesh information                      " <<endl
       << "##########################################################"<<endl<<endl;
        cout << "   minimum mesh size of fluid: hmin = " << min_mesh_size << endl;
        cout << "   maximum mesh size of fliud: hmax = " << max_mesh_size << endl;
 	    }
      else{
        cout <<"**************************************warning***********************************"<<endl;
        cout << "domain_mesh_file has to be specified" << endl;
        cout <<"********************************************************************************"<<endl;        
        libmesh_error();
      } 
  }
  //print mesh info
  if(input_file("print_info", false)){
  	mesh.print_info();
  	cout << endl << endl;
  }
  //===============> (2/4) Create periodic box 
  cout<<"==>(2/4) Create periodic box" <<endl;
  const Point bbox_pmin(XYZAB[0], XYZAB[2], XYZAB[4]);
  const Point bbox_pmax(XYZAB[1], XYZAB[3], XYZAB[5]);
  // construct PMPeriodicBoundary class using info above
  PMPeriodicBoundary pm_periodic_boundary(bbox_pmin, bbox_pmax, periodicity, inlet, inlet_pressure);

  //===============> (3/4) Create polymer chain class
  cout<<"==>(3/4) Create polymer_chain object";
  // Read the polymer data from the local input file
  const unsigned int chain_id = 0;
  PolymerChain polymer_chain(chain_id, pm_periodic_boundary);
  std::ostringstream pfilename;
  if(restart)
  {
	pfilename << point_particle_model<<"_data_restart_"<< restart_step << ".vtk";
	cout <<"-------------> read "<<point_particle_model<<" data from "<<pfilename.str()<< "in restart mode"<<endl;
    		polymer_chain.read_data_vtk(pfilename.str());
  } 
  else
  {
	//generate polymer chains using COPSS 
 		if (point_particle_model=="bead") pfilename << "bead_data_"<<Nb<<"_beads.in";
		else if (point_particle_model=="polymer_chain") pfilename << "polymer_chain_data_"<<nChains<<"_chains.in";
 		cout<<endl<<"--------------> skip generating datafile, will read in existed pizza file: "<<pfilename.str()<<endl;
  	polymer_chain.read_data_pizza(pfilename.str());
		cout 	  <<"--------------> Polymer_chain class is built!" << endl;
		comm_in.barrier();
  }
  pfilename.str(""); pfilename.clear();
  comm_in.barrier();
// Construct PointMesh object from the polymer
  cout <<"==>(4/4) Create point_mesh object, including mesh, polymer_chain and search radius" <<endl;
  const Real search_radius_p = 4.0/alpha;
  const Real search_radius_e = 0.5*max_mesh_size + search_radius_p;
  PointMesh<3> point_mesh(mesh, polymer_chain, search_radius_p, search_radius_e);
  point_mesh.add_periodic_boundary(pm_periodic_boundary);
  point_mesh.reinit();
  cout<<"--------------> Reinit point mesh object, finished!"<<endl;

  if(input_file("print_info",false)){
  	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
  	cout << "### The point-mesh info:\n";
  	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
  	cout << "Total number of point particles: " << point_mesh.num_particles() << "\n";
  	cout <<"search_radius_p = "<<search_radius_p <<", search_radius_e = "<<search_radius_e <<endl<<endl<<endl;
  	if(comm_in.rank()==0){
    		point_mesh.print_point_info();
  	}
  }

  
 /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create an equation systems and ParticleMeshLinearImplicitSystem "Stokes",
   and add variables: velocity (u, v, w) and pressure p. To satisfy LBB condition,
   (u, v, w):  second-order approximation
   p :         first-order basis
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  cout<<endl<<"============================3. Create an equation_systems object(of type 'EquationSystems')============================"<<endl<<endl;
  cout << "==>(1/8) Initialize equation_systems object using the 'mesh' we created before" <<endl;
  EquationSystems equation_systems (mesh);
  
  cout << "==>(2/8) Add 'Stokes' system (of PMLinearImplicitSystem) to the 'equation_systems'" <<endl;
  PMLinearImplicitSystem& system = equation_systems.add_system<PMLinearImplicitSystem> ("Stokes");
  cout << "==>(3/8) Add variables to 'Stokes' system" <<endl;
  unsigned int u_var = 0, v_var = 0, w_var = 0;
  u_var = system.add_variable ("u", SECOND);
  v_var = system.add_variable ("v", SECOND);
  if(dim==3)  w_var  = system.add_variable ("w", SECOND);
  const unsigned int p_var = system.add_variable ("p", FIRST);
  cout<<"--------------> Add variable 'u', 'v', 'w', 'p', finished!"<<endl;
  // attach the particle-mesh system 
  cout << "==>(4/8) Attach point_mesh to the system" <<endl;
  system.attach_point_mesh(&point_mesh);


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Add periodic boundary conditions for the system.
   For the side number of a box, refer to libMesh::Hex27::build_side()
   (https://)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  cout<<"==>(5/8) Add period boundary conditions to 'Stokes' system"<<endl;
  cout << "--------------> Get dof_map of 'Stokes' system"<<endl;
  DofMap& dof_map = system.get_dof_map();
  cout << "--------------> Add periodicBoundary object to 'dof_map'"<<endl;

  /*** set PBC in x-direction ***/
  if (periodicity[0])
  {
    PeriodicBoundary pbcx(RealVectorValue(XYZAB[1]-XYZAB[0], 0., 0.));
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
    
    cout<<"--------------> Set PBC in x direction (for 'u','v','w') finished!"<<endl;
    // check
    if (search_radius_p>=(XYZAB[1]-XYZAB[0])/2.)
    {
      cout<<endl<<"****************************** warning: ********************************"<<endl;
      cout<<"**** The search radius is larger than the domain length in x direction!"<<endl;
      cout<<"**** search radius = "<<search_radius_p<<", domain size Lx ="<< (XYZAB[1]-XYZAB[0])/2.<<endl;
      cout<<"************************************************************************"<<endl<<endl<<endl;
    }
  }
  /*** set PBC in y-direction ***/
  if (periodicity[1])
  {
    PeriodicBoundary pbcy(RealVectorValue(0., XYZAB[3]-XYZAB[2], 0.));
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
    
    cout<<"--------------> Set PBC in y direction (for 'u','v','w'), finished!"<<endl;
    // check
    if (search_radius_p>=(XYZAB[3]-XYZAB[2])/2.)
    {
      cout<<endl<<"****************************** warning: ********************************"<<endl;
      cout<<"**** The search radius is larger than the domain length in y direction!"<<endl;
      cout<<"**** search radius = "<<search_radius_p<<", domain size Ly ="<< (XYZAB[3]-XYZAB[2])/2.<<endl;
      cout<<"************************************************************************"<<endl<<endl<<endl;
    }

  }
  
  
  /*** set PBC in z-direction ***/
  if (periodicity[2])
  {
    PeriodicBoundary pbcz(RealVectorValue(0., 0., XYZAB[5]-XYZAB[4]));
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
    
    cout<<"--------------> Set PBC in z direction (for 'u','v','w'), finished!"<<endl;
    // check
    if (search_radius_p>=(XYZAB[5]-XYZAB[4])/2.)
    {
      cout<<endl<<"****************************** warning: ********************************"<<endl;
      cout<<"**** The search radius is larger than the domain length in z direction!"<<endl;
      cout<<"**** search radius = "<<search_radius_p<<", domain size Lz ="<< (XYZAB[5]-XYZAB[4])/2.<<endl;
      cout<<"************************************************************************"<<endl<<endl<<endl;
    }
  }

    
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize the Preconditioning matrix for saddle point problems if required.
   Initialize the equation system and zero the preconditioning matrix
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if( user_defined_pc ){
  	 cout<<"==> (If user_defined_pc) Add preconditioner matrix to 'Stokes' system"<<endl;
	 system.add_matrix("Preconditioner");
  }	
  /* Initialize the data structures for the equation system. */
  cout<<"==>(6/8) Init equation_systems (libmesh function, to init all systems in equation_systems)"<<endl;
  equation_systems.init ();
  
  // zero the PC matrix, which MUST be done after es.init()
  if( user_defined_pc ) {
  	cout<<"==> (If user_defined_pc) Zero preconditioner matrix"<<endl;
	system.get_matrix("Preconditioner").zero();
  }
  cout<<"--------------> Equation systems are initialized:\n"<<std::endl;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   set parameters of equations systems
   * NOTE: some of these parameters are not used if we use non-dimensional formulation
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  cout<<"==>(7/8) Set parameters of equation_systems"<<endl;
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = max_linear_iterations;
  equation_systems.parameters.set<Real> ("linear solver rtol") = linear_solver_rtol;
  equation_systems.parameters.set<Real> ("linear solver atol") = linear_solver_atol;
  equation_systems.parameters.set<bool>    ("user_defined_pc") = user_defined_pc;
  equation_systems.parameters.set<bool>     ("schur_user_ksp") = schur_user_ksp;
  equation_systems.parameters.set<Real>("schur_user_ksp_rtol") = schur_user_ksp_rtol;
  equation_systems.parameters.set<Real>("schur_user_ksp_atol") = schur_user_ksp_atol;
  equation_systems.parameters.set<string>    ("schur_pc_type") = schur_pc_type;
  equation_systems.parameters.set<StokesSolverType> ("solver_type") = solver_type;
  equation_systems.parameters.set<Real>              ("alpha") = alpha;
  equation_systems.parameters.set<Real>         ("kBT")        = kBT;
  equation_systems.parameters.set<Real>   ("fluid mesh size")  = min_mesh_size;
  equation_systems.parameters.set<Real>       ("viscosity_0")  = muc;
  equation_systems.parameters.set<Real>               ("br0")  = 1.0;
  equation_systems.parameters.set<Real>               ("bk")   = bk;
  equation_systems.parameters.set<Real>               ("q0")   = q0;
  equation_systems.parameters.set<Real>       ("bead radius")  = Rb;
  equation_systems.parameters.set<Real>              ("drag")  = drag_c;
  equation_systems.parameters.set<Real>                ("tc")  = tc;
  equation_systems.parameters.set<Real>               ("Nks")  = Nks;
  equation_systems.parameters.set<Real>               ("Ss2")  = Ss2;
  equation_systems.parameters.set<string> ("particle_type")  = particle_type;
  equation_systems.parameters.set<string> ("point_particle_model") = point_particle_model;
  // Attach force fields
  equation_systems.parameters.set<std::vector<string>> ("pp_force_types") = pp_force_types;
  for (int i=0; i<num_pp_force; i++) equation_systems.parameters.set<std::vector<Real>> (pp_force[i].first) = pp_force[i].second;
  equation_systems.parameters.set<std::vector<string>> ("pw_force_types") = pw_force_types;
  for (int i=0; i<num_pw_force; i++) equation_systems.parameters.set<std::vector<Real>> (pw_force[i].first) = pw_force[i].second;
  equation_systems.parameters.set<string> ("test_name") = test_name;
  equation_systems.parameters.set<string> ("wall_type") = wall_type;
  equation_systems.parameters.set<std::vector<Real>> (wall_type) = wall_params;
 
 
  /* Print information about the mesh and system to the screen. */
  //mesh.print_info();
  if(input_file("print_info",false)){
	cout << endl <<"--------------> Print equation systems info" <<endl;
  	equation_systems.print_info();
  	cout <<"  System has: "<< mesh.n_elem()<<" elements,\n"
             <<"              "<< mesh.n_nodes()<<" nodes,\n"
             <<"              "<< equation_systems.n_dofs()<<" degrees of freedom.\n"
             <<"              "<< equation_systems.n_active_dofs()<<" active degrees of freedom.\n"
             <<"              "<< point_mesh.num_particles()<<" particles.\n" << std::endl;
  }//end if (print_info)

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Init the force field and attach it to the PMLinearImplicitSystem pm_system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  cout<<"==>(8/8) Attach force_field to 'stokes' system"<<endl;
  ForceField force_field(system);
  system.attach_force_field(&force_field);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Start time marching and evolve dynamics of particles
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout<<endl<<"============================4. Start moving particles ============================"<<endl<<endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Parameters for dynamic process
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const unsigned int NP     = point_mesh.num_particles();
  const unsigned int n_vec  = dim*NP;
  bool  cheb_converge;
  Real eig_min = 0., eig_max = 0., real_time = restart_time;

  const Real  max_spring_len = q0/Rb;     // non-dimensional max spring length
  const Real  hmin = equation_systems.parameters.get<Real>("fluid mesh size");

  // Get a better conformation of polymer chains before simulation.
   if(point_particle_model == "polymer_chain")
   {
     bool chain_broken = polymer_chain.check_chain(max_spring_len);
     if(chain_broken) {
         cout << "   ********** warning: Polymer chain is broken in the initial data input" << endl;
         cout << "   ********** warning: bead position is corrected by scaling the chain length and moving the particle according to periodicity" << endl;      
     }
   }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute undisturbed velocity field without particles.
   NOTE: We MUST re-init particle-mesh before solving Stokes
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout<<"==>(1/3) Compute the undisturbed velocity field"<<endl;
  system.reinit_system();
  bool reinit_stokes = true;
  system.solve_stokes("undisturbed",reinit_stokes);
  UniquePtr<NumericVector<Real>> v0_ptr = system.solution->clone(); // backup v0
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   write out the equation systems if write_es = true at Step 0 (undisturbed field)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const std::string out_filename   = "output_pm_system";
  ExodusII_IO*  exodus_ptr = new ExodusII_IO(mesh);
  if(write_es && restart==false)
  {
    //system.add_local_solution(); // Don't add local solution for undisturbed system!
#ifdef LIBMESH_HAVE_EXODUS_API
    exodus_ptr->write_equation_systems(out_filename+".e", equation_systems);
#endif
  }
  
  /* output particle data at the 0-th step in the VTK format */
  if(restart==false)
  {
    if(comm_in.rank()==0){
	     if(point_particle_model == "polymer_chain") polymer_chain.write_polymer_chain("output_polymer_0.vtk");
	     else polymer_chain.write_bead("output_bead_0.csv");
    }
  }

  cout<<"==>(2/3) Prepare RIN & ROUT in binary format at step 0"<<endl;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create vectors and Shell Mat for use:
   U0:          particle velocity vector;
   R0/R_mid:    particle position vector;
   dw/dw_mid:   random vector;
   RIN/ROUT:    the initial and intermediate particle postion vector for msd output
   RIN will not change, and ROUT excludes pbc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec             U0, R0, R_mid, RIN,ROUT, dw, dw_mid;
  Mat             M;
  PetscRandom     rand_ctx;
  PetscViewer     viewer;
  PetscScalar     coef = 0.0;
  BrownianSystem brownian_sys (equation_systems);
  brownian_sys.init_petsc_random(&rand_ctx);
  brownian_sys._create_shell_mat(n_vec, &M);
  brownian_sys._create_petsc_vec(n_vec,&R0);
  VecDuplicate(R0,&U0);
  VecDuplicate(R0,&R_mid);
  VecDuplicate(R0,&dw_mid);
  brownian_sys.extract_particle_vector(&ROUT,"coordinate","extract");
  VecDuplicate(ROUT,&RIN);
  VecCopy(ROUT,RIN);  // RIN = ROUT = the initial position vector
  brownian_sys.set_std_random_seed(random_seed);  // random seed
 
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


  /* Calculate the intial center of mass */
  const std::vector<Point> center0 = brownian_sys.center_of_mass(RIN);

  if(!restart)
  {
    /* Output mean square displacement and radius of gyration at step 0 (the origin) */
    brownian_sys.output_statistics_step0(out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag, RIN);
  }

  cout<<"==>(3/3) Start calculating dynamics and advancing time steps"<<endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Advancing in time. Fixman Mid-Point algorithm
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const unsigned int istart = restart_step;
  const unsigned int iend   = restart_step + nstep;
  unsigned int o_step = restart_step;  // output step
  std::vector<Real> vel0(n_vec,0.);
  std::vector<Real> vel1(n_vec,0.);
  for(unsigned int i=istart; i<iend; ++i)
  {
    cout << endl<<"   Starting Fixman Mid-Point algorithm at step " << i+1 << "..."<<endl;
 
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
    brownian_sys.vector_transform(vel1, &U0, "forward");
 
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
        exodus_ptr->write_timestep(out_filename+".e",system.get_equation_systems(),o_step,o_step);
        #endif
      } // end if (write es)   
    } // end if (i % write_interval == 0 )
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adaptive time step.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Real dt = 0;
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
        cout << "       ##############################################################################################################" << endl
             << "       # Max velocity magnitude is " << vp_max << endl
             << "       # Min velocity magnitude is " << vp_min << endl
             << "       # minimum fluid mesh size = " << hmin << endl
             << "       # The adaptive time increment at step "<< i+1 << " is dt = " << dt<<endl
             << "       # (with Brownian) adapting_time_step = max_dr_coeff * bead_radius / (max_bead_velocity at t_i)" << endl
             << "       # (without Brownian) adapting_time_step = max_dr_coeff * fluid_mesh_size_min / (max_bead_velocity at t_i) "<<endl      
             << "       ##############################################################################################################" << endl;
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
        cout << "       ##############################################################################################################" << endl
             << "       # The fixed time increment at step "<< i+1 << " is dt = " << dt<<endl
             << "       # (With Brownian) fixed_time_step = max_dr_coeff * bead_radius / 1.0"<<endl
             << "       # (Without Brownian) fixed_time_step = max_dr_coeff * fluid_mesh_size_min / 1.0 "<<endl
             << "       ##############################################################################################################" << endl;
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
      brownian_sys.std_random_vector(0.0,std_dev,"gaussian",&dw);
      brownian_sys._vector_mean_variance(dw, mean_dw, variance_dw);
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
          brownian_sys.compute_eigenvalues(eig_min,eig_max,tol_eigen);
        }
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the Brownian displacement B^-1 * dw using Chebyshev approximation.
         Here dw is both input and output variables, so it will be changed.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        //PetscPrintf(PETSC_COMM_WORLD,
        //            "--->eig_min = %f, eig_max = %f, tol_cheb = %f, max_n_cheb = %d\n",
        //            eig_min,eig_max,tol_cheb,max_n_cheb);
        cheb_converge = brownian_sys.chebyshev_polynomial_approximation(max_n_cheb,
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
        cout << "****** Warning: Chebysheve failed to converge at step " <<i+1<<endl;
      }

      // Magnify the spectral range by a factor (1.05 by default).
      eig_max *= eig_factor; eig_min /= eig_factor;
      
      
      // Compute dw_mid = D*B^-1*dw, which can be obtained by solving the Stokes
      brownian_sys.hi_ewald(M,dw,dw_mid);  // dw_mid = D * dw

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Particle coordinate vector R0.
       Move the particle R_mid = R0 + 0.5*(U0+U1)*dt (deterministic)
       and R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw     (stochastic)
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      brownian_sys.extract_particle_vector(&R0,"coordinate","extract");
      VecWAXPY(R_mid,0.5*dt,U0,R0);  // R_mid = R0 + 0.5*dt*(U0+U1)  (R0 and U0 do NOT change)
      coef = 0.5;                    // coefficient. sqrt(2) is introduced when generating dw
      VecAXPY(R_mid,coef,dw_mid);    // R_mid = R_mid + 0.5*sqrt(2)*D*B^-1*dw
      brownian_sys.extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords
   
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Check and correct beads' position at the midpoint
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      force_field.check_walls(); // check pbc and inpenetrable wall
      if(point_particle_model == "polymer_chain")
      {
        bool chain_broken = polymer_chain.check_chain(max_spring_len);
        if(chain_broken) {
          cout << "   ********** warning: Polymer chain is broken in midpoint at the step " <<i+1<<endl;
          cout << "   ********** warning: bead position is corrected by scaling the chain length and moving the particle according to periodicity" << endl;      }
      }
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
      brownian_sys.vector_transform(vel1, &U0, "forward"); // (U0+U1)_mid
      brownian_sys.hi_ewald(M,dw,dw_mid);  // dw_mid = D_mid*dw, where dw=B^-1*dw computed above
   
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       the mid-point to the NEW point, and update the particle coordinates
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      //cout << "Update from mid-point to the NEW particle coordinates at step " <<i+1<<endl;
      VecWAXPY(R_mid,dt,U0,R0);         // R_mid = R0 + dt*U0_mid
      VecAXPY(R_mid,2.0*coef,dw_mid); // R_mid = R_mid + sqrt(2)*D_mid*B^-1*dw
      brownian_sys.extract_particle_vector(&R_mid,"coordinate","assign");
  
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Check and correct the beads' position again after the midpoint update
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      force_field.check_walls(); // check pbc and inpenetrable wall
      if(point_particle_model == "polymer_chain")
      {
        bool chain_broken = polymer_chain.check_chain(max_spring_len);
        if(chain_broken) {
          cout << "   ********** warning: Polymer chain is broken at the step " <<i+1<<endl;
        }
      }

      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT,dt,U0);            // ROUT = ROUT + dt*U0_mid
      VecAXPY(ROUT,2.0*coef,dw_mid);  // ROUT = ROUT + sqrt(2)*D_mid*B^-1*dw
    } // end if Brownian
    else{ // if without Brownian

      // Move the particle R_mid = R0 + (U0+U1)*dt (deterministic)
      brownian_sys.extract_particle_vector(&R0,"coordinate","extract");
      VecWAXPY(R_mid,dt,U0,R0);  // R_mid = R0 + dt*Utotal (U0 is actually Utotal)
      brownian_sys.extract_particle_vector(&R_mid,"coordinate","assign"); // Update mid-point coords
 
      // Check and correct beads' position at the midpoint
      force_field.check_walls(); // check pbc and inpenetrable wall
      if(point_particle_model == "polymer_chain")
      {
        bool chain_broken = polymer_chain.check_chain(max_spring_len);
        if(chain_broken) {
          cout << "   ********** warning: Polymer chain is broken in midpoint at the step " <<i+1<<endl;
          cout << "   ********** warning: bead position is corrected by scaling the chain length and moving the particle according to periodicity" << endl;      }
      }

      // Update ROUT (position vector excluding pbc) at the i-th step
      VecAXPY(ROUT,dt,U0); // ROUT = ROUT + dt*U0_mid
 
    } // end else (without_brownian)

    // Update the time
    real_time += dt;
    
    if(i%write_interval==0){
      o_step++; 

      if(comm_in.rank()==0){

        /*---------------------------------------------------------------------------------------
         * output polymer chain / bead data in the VTK format at step i
        -----------------------------------------------------------------------------------------*/
      	if(point_particle_model == "polymer_chain"){
                 oss << "output_polymer_" << o_step << ".vtk";
      	   polymer_chain.write_polymer_chain( oss.str() );
      	}
      	else{
      	   oss << "output_bead_" << o_step <<".csv";
      	   polymer_chain.write_bead(oss.str());
      	} // end else

        /*----------------------------------------------------------------------------------------------------
         * Output mean square displacement, radius of gyration, chain stretch, and center of mass at step i
         --------------------------------------------------------------------------------------------------- */
        brownian_sys.output_statistics_stepi(out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag,
                                             i, real_time, center0, ROUT);
      } // end if comm_in.rank() == 0

      oss.str(""); oss.clear();
      /*----------------------------------------------------------------------------------------------------
       * Write out ROUT for restart mode at step i
      ----------------------------------------------------------------------------------------------------*/
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vector_ROUT.dat",FILE_MODE_WRITE,&viewer);
      VecView(ROUT,viewer);

    }
  } //end this time step

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Destroy and return
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
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

  time(&rawtime);
  timeinfo = localtime ( &rawtime );
  cout <<endl
       <<"----------------------------------------------------------------------"<<endl
       <<"The current date/time is: "<<asctime(timeinfo)
       <<"The simulation is finished and Destroy all the PETSc objects."
       <<"----------------------------------------------------------------------"<<endl<<endl;  



  // return
  return 0;
}
