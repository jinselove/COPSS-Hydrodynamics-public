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

  void read_data(std::string control_file);


  /*
   * Read test_name
   */
  virtual void read_test_name();

  /*
   * Read physical parameters 
   */
  virtual void read_physical_parameters();

  /*
   * Read Geometry infomation
   */

  virtual void read_geometry();

  /*
   * Read mesh
   */
  virtual void read_mesh();

  /*
   * read force types
   */
  virtual void read_force();

  /*
   * read GGEM info
   */
  virtual void read_ggem();

  /*
   * read Stokes Solver  
   */
  virtual void read_stokes_solver();

  /*
   * read Chebyshev info
   */
  virtual void read_chebyshev();

  /*
   * read run time info 
   */
  virtual void read_run();


  //protected variables

  // PETSC MPI communicator
  Parallel::Communicator _comm_in;
  // control file object
  GetPot _input_file;
  // test name
  std::string _test_name;
  // physical parameters
  std::string _particle_type;
  Real _viscosity, _drag_c, _Rb, _Db;

};

} // end namespace libMesh
