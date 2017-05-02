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


#pragma once

#include <stdio.h>

#include "libmesh/reference_counted_object.h"
#include "libmesh/parallel_object.h"

#include "polymer_chain.h"
#include "rigid_particle.h"



namespace libMesh
{
  
  /*
   * This class defines a particle-chain object, in which particles
   * are respresented by finite element mesh, and the chains are described 
   * by bead-spring model
   */
  
  
class ParticleChain :  public ReferenceCountedObject<ParticleChain>,
                       public ParallelObject
{
public:
  
  // Constructor
  ParticleChain(const std::size_t n_particles,
                const std::size_t n_chains,
                const Parallel::Communicator &comm_in);
  
  
  // Destructor
  ~ParticleChain();
  
  
  /*
   * Read the data of chain from local file, whose data structure is:
   * Nb
   * # type x0 y0 z0 a0 b0 c0 theta0
   * # type x1 y1 z1 a1 b1 c1 theta1
   * ......
   *
   * particle_type_id: type id to identify the point is a particle
   */
  void read_data(const std::string& filename,
                 const std::string& vmesh_file,
                 const std::string& smesh_file,
                 const std::string& mesh_type,
                 const std::size_t  particle_type_id = 1);
  
  
  /*
   * read partilces
   */
  void read_particles(const std::string& filename,
                      const std::string& vmesh_file,
                      const std::string& smesh_file,
                      const std::string& mesh_type,
                      const std::size_t  particle_type_id = 1);
  
  
  /*
   * read chains
   */
  void read_chains(const std::string& filename,
                   const std::size_t  particle_type_id = 1);
  
  
private:
  // number of particles
  std::size_t _n_particles;
  
  
  // number of particles
  std::size_t _n_chains;
  
  /*
   * If we use std::vector<PolymerChain> or std::vector<RigidParticle>,
   * vector's resize() function doesn't know how much memory to allocate
   * for the objects. Therefore, we use pointers instead.
   */
  
  // particles
  std::vector<RigidParticle*> _particles;

  // polymer chains (can be more than one chain)
  std::vector<PolymerChain*> _chains;
  
};  // end of class
  
}  // end of namespace
