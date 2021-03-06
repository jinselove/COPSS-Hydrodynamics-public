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
#include "libmesh/perf_log.h"

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
  bool print_info;
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
  unsigned int num_pp_force;
  std::vector<std::string> pp_force_type;
  std::vector<ForceField::type_force> pp_force;
  unsigned int num_pw_force;
  std::vector<std::string> pw_force_type;
  std::vector<ForceField::type_force> pw_force;

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
  Real max_dr_coeff = 0.1; // max displacement per step
  bool restart; // if restart
  std::size_t restart_step; // restart step
  Real restart_time; // real time at restart step
  unsigned int nstep; // totol number of steps to run
  unsigned int write_interval; // output file write interval
  bool write_es, out_msd_flag, out_stretch_flag, out_gyration_flag, out_com_flag;
  std::ostringstream oss;

  // mesh
  SerialMesh* mesh;
  PointMesh<3>* point_mesh;
  //std::unique_ptr<SerialMesh> mesh;
  Real min_mesh_size, max_mesh_size;
  Real search_radius_p, search_radius_e;

  //periodic boundary
  PMPeriodicBoundary* pm_periodic_boundary;
  //std::unique_ptr<PMPeriodicBoundary> pm_periodic_boundary;

  //force field
  ForceField* force_field;
  // equation system
  unsigned int u_var, v_var, w_var, p_var;

  //integrate
  // paramters for dynamic process;
  bool reinit_stokes;
  unsigned int NP;
  unsigned int n_vec;
  Real hmin;
  bool cheb_converge;
  Real eig_min = 0, eig_max = 0;
  Real real_time;
  std::string out_system_filename = "output_pm_system";
  UniquePtr<NumericVector<Real>> v0_ptr;
  ExodusII_IO* exodus_ptr;

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
  BrownianSystem* brownian_sys; 
  std::vector<Point> center0;

  //variables in integrator
  unsigned int o_step;
  std::vector<Real> vel0;
  std::vector<Real> vel1;

  





  /*!
   * class constructor
   */
  Copss(const CopssInit& init);

  /*!
   * class destructor
   */
  ~Copss();
  
  /*! Check libmesh support
   *
   * PETSC or laspack
   * SLEPC
   * AMR
   * at least 2D
  */
  int check_libmesh();

  /*!
   * Print out start time
   */ 
  void start_time(struct tm * timeinfo);

  /*!
   * Print out end time
   */ 
  void end_time(struct tm * timeinfo);

  /*!
   * Build and initialized an Equation system
   * step 1: read_input()
   * step 2: create_object_mesh() 
   * step 3: create_equation_systems()
   */ 
  EquationSystems init_system(std::string input_file); // including the following 10 steps

  /*!
   * Read input and create necessary objects
   * step 1: read_system_info()
   * step 2: read_physical_info()
   * step 3: read_particle_info()
   * step 4: read_domain_info()
   * step 5: read_force_info()
   * step 6: read_ggem_info()
   * step 7: read_stokes_solver_info()
   * step 8: read_chebyshev_info()
   * step 9: read_run_info()
   */
  void read_input();

  
  /*! 
   * Create object_mesh 
   * this object_mesh can be point_mesh or particle_mesh
   * step 1: create_domain_mesh()
   * step 2: create_periodid_boundary()
   * step 3: create_object()
   * step 4: create object_mesh by combining domain mesh, object and periodic boundary
   */
  virtual void create_object_mesh() = 0;

  /*! 
   * Create equation_system object
   * step 1: initialize equation_systems using mesh;
             add a PMLinearImplicitSystem, system, to equation_systems;
             add variables
   * step 3: attach_object_mesh()
   * step 4: attach_period_boundary()
   * step 5: initialize equation systems and preconditioning matrix
   * step 5: set_parameters() for equation_systems
   * step 6: Initialize the force field
   */
  EquationSystems create_equation_systems();

  /*!
   * integrater
   * step 1: 
  */
  virtual void run(EquationSystems& equation_systems) = 0;

  /*!
   * destroy PETSC objects
   */
  void destroy(); 


protected:
  /*!
   * Steps for read_input() 
   */

  void read_system_info(); 
  void read_physical_info(); 
  virtual void read_particle_info() = 0;
  void read_domain_info(); 
  void read_force_info(); 
  void read_ggem_info(); 
  void read_stokes_solver_info(); 
  void read_chebyshev_info(); 
  void read_run_info(); 

  /*!
   * Steps for create_object_mesh()
   */
  void create_domain_mesh();
  void create_periodic_boundary();
  virtual void create_object() = 0;

  /*!
   * Steps for create_equation_systems
   */
  virtual void attach_object_mesh(PMLinearImplicitSystem& system) = 0;
  void attach_period_boundary(PMLinearImplicitSystem& system);
  virtual void set_parameters(EquationSystems& equation_systems) = 0;

  /*!
   * Steps for integrate()
   */
  /*!
   * Compute undisturbed velocity field without particles.
   * NOTE: We MUST re-init particle-mesh before solving Stokes
   */ 
  void solve_undisturbed_system(EquationSystems& equation_systems);

  /*!
  * Create vectors and Shell Mat for use
  * Create brownian system
  */ 
  void create_brownian_system(EquationSystems& equation_systems);

  /*!
   * Integrate particle motions using Fixman's midpoint scheme
   */
  void fixman_integrate(EquationSystems& equation_systems, unsigned int i);
  // update object positions due to PBS
  virtual void update_object(std::string stage) = 0;

};

} // end namespace libMesh
