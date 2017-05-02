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


//
//  stokes_solver.h
//  
//
//  Created by Xujun Zhao on 2/17/16.
//
//

#pragma once

#include "libmesh/petsc_macro.h"  // ahead of LIBMESH_HAVE_PETSC
#include "libmesh/libmesh_config.h"


// make sure libMesh has been complied with PETSc
#ifdef LIBMESH_HAVE_PETSC          // 11111111111111111111111111111111


// C++ includes
#include <iostream>
//#include <memory>

// Local includes libmesh headers
#include "libmesh/equation_systems.h"
#include "libmesh/parallel_object.h"
#include "libmesh/reference_counted_object.h"


// include PETSc KSP solver
EXTERN_C_FOR_PETSC_BEGIN
#  include <petscksp.h>
EXTERN_C_FOR_PETSC_END


using libMesh::EquationSystems;
using libMesh::ReferenceCountedObject;
using libMesh::ParallelObject;
using libMesh::Real;



/* this class performs schur complement reduction type solve
 for the saddle point problems arised from mixed finite
 element methods
 */



enum StokesSolverType
{
  superLU_dist,
  field_split,
  user_define
};



class StokesSolver : public ReferenceCountedObject<StokesSolver>,
                     public ParallelObject
//class SchurComplementSolver
{
public:
  // constructor
  StokesSolver(EquationSystems& es_stokes);
               
               
  // constructor
  StokesSolver(EquationSystems& es_stokes,
               const StokesSolverType solver_type);
  
  
  // Destructor
  ~StokesSolver();
  
  
  /*
   * Init the KSP solver:
   * The system matrix needs to be assembled before calling this init function! 
   */
  void init_ksp_solver();
  
  
  /*
   * Return if the ksp solver is initialized or not.
   */
  const bool is_ksp_initialized() const { return _is_init;  }
  
  
  /*
   * Set the solver type
   */
  void set_solver_type(const StokesSolverType solver_type);
  
  
  
  /*
   * solve the equation system Ax = b
   */
  void solve();
  
  
  /*
   * Build IS for u/p for PC field split
   *    -called by: solve()
   */
  void build_is(IS *is_v,
                IS *is_p);
  
  
  /*
   * Set up the PC for the schur complement reduction algorithm
   *    -called by solve()
   */
  void setup_schur_pc(KSP ksp,
                      IS is_v,
                      IS is_p,
                      Mat *pmat,
                      const bool userPC,
                      const bool userKSP);
  
  
  /*                                ^
   * Schur complement approximation S used as a PC for S*y1=y2, for example:
   *     S1 = A11 - A10 diag(A00)^(-1) A01;
   *     Sa = 1/v * Mp (pressure matrix, v is kinematic viscosity)
   *    -called by solve()
   */
  void setup_approx_schur_matrix(IS is_v,
                                 IS is_p,
                                 Mat *pmat);
  
  
  
  /*
   * Return the number of linear iterations required to solve Ax=b
   */
  const unsigned int n_linear_iterations() const;
  
  
  /*
   * Return the final residual of the linear solve
   */
  const Real final_linear_residual() const;
  
  
  /*
   * Petsc View
   */
  void petsc_view_is(IS is_p) const;
  void petsc_view_vector(Vec vector) const;
  void petsc_view_matrix(Mat matix) const;
  
  
private:
  
  // EquationSystems* _es;
  EquationSystems& _equation_systems;
  
  // the type of the solver used for Stokes
  StokesSolverType _solver_type;
  
  // solver relative tolerance
  PetscReal _rtol;
  
  // solver absolute tolerance
  PetscReal _atol;
  
  // (iterative) solver maximum iteration
  PetscInt _max_it;

  // KSP sover
  KSP _ksp;
  
  // IS pointers for velocity and pressure, respectively
  IS _is_v;
  IS _is_p;
  
  // preconditioning matrix for Schur Complement
  Mat _schur_pmat;
  
  // Label if the system is initialized.
  // If not, the destructor cannot destroy PETSc objects.
  bool _is_init;
};  // end of class SchurComplementSolver
#endif // #ifdef LIBMESH_HAVE_PETSC // 11111111111111111111111111111111


