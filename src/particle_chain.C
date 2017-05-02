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



#include <fstream>

#include "libmesh/libmesh_common.h"

#include "pm_toolbox.h"
#include "particle_chain.h"




namespace libMesh
{




// ===========================================================================
ParticleChain::ParticleChain(const std::size_t n_particles,
                             const std::size_t n_chains,
                             const Parallel::Communicator &comm_in)
: ParallelObject(comm_in),
  _n_particles(n_particles), _n_chains(n_chains)
{
  // init the containers
  _particles.resize(n_particles);
  _chains.resize(n_chains);
}




// ===========================================================================
ParticleChain::~ParticleChain()
{
  // Free the memory for RigidParticle
  for(std::size_t i=0; i<_particles.size(); ++i){
    delete _particles[i];
  }
  
  // FIXME: Do we need to free the memory for PolymerChain(PointParticles)?
}




// ===========================================================================
void ParticleChain::read_data(const std::string& filename,
                              const std::string& vmesh_file,
                              const std::string& smesh_file,
                              const std::string& mesh_type,
                              const std::size_t  particle_type_id)
{
  START_LOG ("read_data()", "ParticleChain");
  
  std::cout <<"\n### Read particle-chain filename = "<<filename <<std::endl;
  this->read_particles(filename,vmesh_file,smesh_file,mesh_type, particle_type_id);
  
  std::cout << "Reading polymer chain data from "<<filename<<" is completed!\n\n";
  
  STOP_LOG ("read_data()", "ParticleChain");
}



// ===========================================================================
void ParticleChain::read_particles(const std::string& filename,
                                   const std::string& vmesh_file,
                                   const std::string& smesh_file,
                                   const std::string& mesh_type,
                                   const std::size_t  particle_type_id)
{
  START_LOG ("read_particles()", "ParticleChain");
  
  /* --------------------------------------------------------------------------
   * Check if the mesh files exist. If the surface mesh file does NOT exist,
   * it will be extracted from the volume mesh file. However, the volume mesh
   * file MUST exist. Otherwise, print the error message.
   * --------------------------------------------------------------------------*/
  const bool smesh_exist = PMToolBox::file_exist(smesh_file);
  const bool vmesh_exist = PMToolBox::file_exist(vmesh_file);
  if( !vmesh_exist )
  {
    printf("***error in read_particles_data(): volume mesh file does NOT exist!");
    libmesh_error();
  }
  
  
  /* --------------------------------------------------------------------------
   * check the existance of the particle-chain input file
   * --------------------------------------------------------------------------*/
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***error in read_particles(): particle-chain input file does NOT exist!");
    libmesh_error();
  }
  
  
  /* --------------------------------------------------------------------------
   * initialize:
   * --------------------------------------------------------------------------*/
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  Real x=0., y=0., z=0.;            // initialize bead coords
  std::size_t bead_id, bead_type;   //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  
  
  /* --------------------------------------------------------------------------
   * read particle data, and generate particle surface mesh for rigid particles
   * --------------------------------------------------------------------------*/
  unsigned int n_points = 0, particle_count = 0;
  infile >>    n_points;  // total number of points (NOT particles!)
  for (std::size_t i=0; i<n_points; ++i)
  {
    /*
     * Read the information of each line
     */
    infile >> bead_id >> bead_type >> x >> y >> z
           >>rot_vec[0] >> rot_vec[1] >>rot_vec[2] >>rot_vec[3] ;
    
    // if this is a particle, construct it!
    if(bead_type == particle_type_id)
    {
      // Build the particle
      Point pt(x,y,z);
      RigidParticle* particle = new RigidParticle(pt, particle_count, this->comm());
      
      // ------- the magnification factor and rotation angles -------
      std::vector<Real> mag_factor(3), angles(3);
      const Real r = 5, h = 10;
      mag_factor[0] = r;  mag_factor[1] = r;  mag_factor[2] = h/2.;
      angles[0]  =  0;  angles[1]  =  90;  angles[2]  =  0;
      // ------------------------------------------------------------
      
      // Read the particle's mesh
      if(mesh_type=="surface_mesh")
      {
        if( (!smesh_exist) || particle_count==0 )
          particle->extract_surface_mesh(vmesh_file,smesh_file);
        
        particle->read_mesh_cylinder(smesh_file, mesh_type, mag_factor, angles);
      }
      else if(mesh_type=="volume_mesh")
      {
        particle->read_mesh_cylinder(smesh_file, mesh_type, mag_factor, angles);
      }
      else
      {
        printf("***error in read_particles_data(): invalid mesh type!");
        libmesh_error();
      }
      
      // Put it to the container
      _particles[particle_count] = particle;
      particle_count++;
      
      // --------------------- test ------------------------------------------
      if(this->comm().rank()==0 )
      {
        printf("ParticleChain::read_particles: x = %f, y = %f, z = %f. \n",x, y, z );
        printf("        rotation vector = (%f, %f, %f). \n",angles[0], angles[1], angles[2] );
        printf("        MPI_rank = %d\n",this->comm().rank() );
        printf("        # of elements is %u\n",particle->mesh().n_elem() );
        printf("        # of nodes is %u\n",particle->mesh().n_nodes() );
      }
      // --------------------- test ------------------------------------------
      
    } // end if(bead_type == particle_type_id)
    
  } // end for i-loop
  
  
  /* --------------------------------------------------------------------------
   * close the file and end the function
   * --------------------------------------------------------------------------*/
  infile.close();
  this->comm().barrier();
  
  STOP_LOG ("read_particles()", "ParticleChain");
}


    
// ===========================================================================
void ParticleChain::read_chains(const std::string& filename,
                                const std::size_t  particle_type_id)
{
  START_LOG ("read_chains()", "ParticleChain");
  
  // Open the local file
  std::ifstream infile;
  infile.open (filename, std::ios_base::in);
  if( !infile.good() )
  {
    printf("***warning: read_data() can NOT read the polymer chain data!");
    libmesh_error();
  }
  
  
  /* --------------------------------------------------------------------------
   * initialize:
   * --------------------------------------------------------------------------*/
  // point_type:  0 - polymer bead point; 1 - tracking point; or user-defined type
  Real x=0., y=0., z=0.;            // initialize bead coords
  std::size_t bead_id0, bead_id, bead_type, chain_id = 0;   //
  std::vector<Real> rot_vec(4); // rotation vector (a,b,c) + theta.
  bool build_new_chain = false;
  
  
  /* --------------------------------------------------------------------------
   * read particle data, and generate particle surface mesh for rigid particles
   * --------------------------------------------------------------------------*/
  unsigned int n_points = 0, particle_count = 0;
  infile >>    n_points;  // total number of points (NOT particles!)
  for (std::size_t i=0; i<n_points; ++i)
  {
    /*
     * Read the information of each line
     */
    infile >> bead_id >> bead_type >> x >> y >> z
           >>rot_vec[0] >> rot_vec[1] >>rot_vec[2] >>rot_vec[3] ;
    
    /* If this is a particle, not a bead on a polymer chain,
       the next point will be the starting bead of the chain */
    if(bead_type == particle_type_id)
    {
      if(bead_id==1) {
        chain_id = 0;
      }
      else{
        chain_id++;
      }
      particle_count++;
    }
    
    // Build the PolymerChain
    PolymerChain* polymer_chain = new PolymerChain(chain_id);
    
    
  }
  
  
  STOP_LOG ("read_chains()", "ParticleChain");
}
  
  
  
  
  
  

} // end of namespace
