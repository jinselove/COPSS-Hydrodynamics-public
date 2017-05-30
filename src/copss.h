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
#pragma once

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
#include "point_particle.h"
#include "particle_mesh.h"
#include "point_mesh.h"
#include "force_field.h"
#include "pm_linear_implicit_system.h"
#include "brownian_system.h"
#include "pm_periodic_boundary.h"
#include "chebyshev.h"
#include "pm_toolbox.h"
#include "polymer_chain.h"
#include "random_generator.h"
#include "stokes_solver.h"
#include "ggem_system.h"
#include "copss_init.h"


/*! This file serves as a template to build\solve a Stokes system
    using our pFE-GgEm code 
*/

// Bring in everything from the libMesh namespace


namespace libMesh{

class Copss
{
public: 
  // PETSC MPI communicator
  Parallel::Communicator comm_in;
  // error message string
  std::string error_msg;
  // output message string
  std::string output_msg;
  // control file object
  GetPot input_file;
  // control file name
  std::string control_fileName;
  // test name
  std::string test_name;
  // physical parameters
  const Real kB = 1.380662E-17;//1.380662E-23(J/K) = 1.3806623E-23(N*m/K) = 1.380662E-17(N*um/K)
  const Real PI = libMesh::pi;
  Real T; // simulation temperature (K)
  Real kBT;// (N*um)
  Real viscosity; // viscosity of the fluid (cp = N*s/um^2)
  Real Rb; // radius of the bead
  Real drag_c; // Drag coefficient (N*s/um)
  Real Db;
  std::string particle_type;

  // characteristic variables
  Real tc; // characteristic time (diffusion time) (s)
  Real uc; // characteristic velocity (um/s)
  Real fc; // characteristic force (N)
  Real muc; // non-dimensional viscosity

  // Geometry information
  unsigned int dim; // dimension of the box
  std::string wall_type; // wall_type (slit or sphere)
  std::vector<Real> wall_params; // (wall parameters; e.g. slit = '-50,50,-50,50,-50,50' or sphere = '10')
  std::vector<bool> periodicity; // periodicity of the box
  std::vector<bool> inlet; // inlet direction of the box
  std::vector<Real> inlet_pressure; // inlet pressure

  // Mesh information
  bool generate_mesh; // flag to generate mesh or load mesh
  std::string domain_mesh_file; // domain mesh filename
  std::vector<unsigned int> n_mesh; // mesh size in all directions

  // Force information
  unsigned int num_pp_forces;
  std::vector<std::string> pp_force_types;
  std::vector<ForceField::type_force> pp_forces;
  unsigned int num_pw_forces;
  std::vector<std::string> pw_force_types;
  std::vector<ForceField::type_force> pw_forces;

  // GGEM information
  Real alpha;

  // Solver information
  int max_linear_iterations;
  Real linear_solver_rtol;
  Real linear_solver_atol;
  bool user_defined_pc;
  bool schur_user_ksp;
  Real schur_user_ksp_rtol;
  Real schur_user_ksp_atol;
  std::string schur_pc_type;
  std::string stokes_solver_type;
  StokesSolverType solver_type; 

  // Chebyshev information
  unsigned int max_n_cheb; // max order of Chebyshev polynomianl 
  Real tol_cheb; // tolerance of chebyshev convergence
  Real eig_factor; // factor to resize eigen value range
  Real tol_eigen; // tolerance of eigenvalue convergence
  bool compute_eigen; // if compute eigen value; initially set it true

  // Run time info
  bool with_brownian; // if consider brownian motion
  std::size_t random_seed;
  bool adaptive_dt; // if use adaptive time step (essential for brownian systems)
  Real dt0; // timestep when brownian system starts (step1)
  Real max_dr_coeff; // max displacement per step
  bool restart; // if restart
  std::size_t restart_step; // restart step
  Real restart_time; // real time at restart step
  unsigned int nstep; // totol number of steps to run
  unsigned int write_interval; // output file write interval
  bool write_es, out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag;

  // mesh
  SerialMesh* mesh;
  Real min_mesh_size, max_mesh_size;

  //periodic boundary
  PMPeriodicBoundary* pm_periodic_boundary;

  // class constructor
  Copss(const CopssInit& init);
  
  /* 
   * Check libmesh support
   * PETSC or laspack
   * SLEPC
   * AMR
   * at least 2D
  */
  int check_libmesh();

  /*
   * Print out start time
   */ 
  void start_time(struct tm * timeinfo);

  /*
   * Print out end time
   */ 
  void end_time(struct tm * timeinfo);

  /*
   * This function will be overriden later in inheritance class
   * This function contains all we need to build a new system
   */ 
  void init_system(std::string input_file); // including the following 10 steps
  // (1/10) read all the input information from "input_file"
  void read_input();
  // (2/10) system name 
  void read_system_info(); 
  // (3/10) physical parameters
  void read_physical_info(); 
  // (4/10) information about particles
  virtual void read_particle_info() = 0;
  // (5/10) domain information, including geometry and mesh 
  void read_domain_info(); 
  // (6/10) force types, including particle-particle and particle-wall
  void read_force_info(); 
  // (7/10) ggem information
  void read_ggem_info(); 
  //(8/10) stokes solver control parameter
  void read_stokes_solver_info(); 
  // (9/10) chebyshev parameter
  void read_chebyshev_info(); 
  // (10/10) running parameter
  void read_run_info(); 


  /*
   * Create object_mesh 
   * this object_mesh can be point_mesh or particle_mesh
   * which be defined in corresponding derived classes.
   */

  // generate or create domain mesh
  void create_domain_mesh();
  // create domain periodic box
  void create_periodic_boundary();
  // create object
  virtual void create_object() = 0;
  // create object_mesh (point_mesh or particle_mesh)
  virtual void create_object_mesh() = 0;




};

} // end namespace libMesh
