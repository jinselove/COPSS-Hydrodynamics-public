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


/*! This file serves as a template to build\solve a Stokes system
    using our pFE-GgEm code 
*/

// Bring in everything from the libMesh namespace


namespace libMesh{

class Copss
{
public: 

  // class constructor 
  Copss(int argc, char** argv);
  
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






protected:

  /*
   * This function will be overriden later in inheritance class
   * This function contains all we need to build a new system
   */ 
  virtual void build_system() = 0;

  /*
   * Read input parameters
   */ 

  //virtual void read_data(std::string control_file) = 0;
  virtual void read_data(std::string control_file) final;

  /*
   * Read test_name
   */
  virtual void read_test_info() final;

  /*
   * Read physical parameters
   */
  virtual void read_physical_info() final;

  /*
   * Read physical parameters (will be overriden in derived class)
   */
  virtual void read_particle_info() = 0;

  /*
   * Read Geometry infomation
   */

  virtual void read_geometry_info() final;

  /*
   * Read mesh
   */
  virtual void read_mesh_info() final;

  /*
   * read force types
   */
  virtual void read_force_info() final;

  /*
   * read GGEM info
   */
  virtual void read_ggem_info() final;

  /*
   * read Stokes Solver  
   */
  virtual void read_stokes_solver_info() final;

  /*
   * read Chebyshev info
   */
  virtual void read_chebyshev_info();

  /*
   * read run time info 
   */
  virtual void read_run_info();


  //protected variables

  // PETSC MPI communicator
  Parallel::Communicator comm_in;
  // error message string
  std::string error_msg;
  // control file object
  GetPot input_file;
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
  std::vector<Real> n_mesh; // mesh size in all directions

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


};

} // end namespace libMesh
