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




#include <iomanip>  // setw | setprecision

#include "libmesh/elem.h"

#include "point_particle.h"


namespace libMesh
{
  


// ======================================================================
PointParticle::PointParticle(const Point pt,
                             const dof_id_type point_id)
: _center(pt), _counter(3,0), _id(point_id),
  _point_type(NOT_DEFINED), _parent_id(-1),
  _processor_id(-1),
  _elem_id(-1),
  _force(3,0.)
{
  // do nothing
}



// ======================================================================
PointParticle::PointParticle(const Point pt,
                             const dof_id_type point_id,
                             const PointType point_type)
: _center(pt), _counter(3,0), _id(point_id),
_parent_id(-1), _point_type(point_type),
_processor_id(-1),
_elem_id(-1),
_force(3,0.)
{
  // do nothing
}


  
// ======================================================================
PointParticle::PointParticle(const Point pt,
                             const dof_id_type point_id,
                             const PointType point_type,
                             const std::vector<Real>& rot_vec)
  : _center(pt), _counter(3,0), _id(point_id),
  _parent_id(-1), _point_type(point_type),
  _processor_id(-1),
  _elem_id(-1),
  _force(3,0.),_orientation(rot_vec)
{
  // do nothing
}
  
  

// ======================================================================
PointParticle::PointParticle(const PointParticle& particle)
{
  _center       = particle._center;
  _counter      = particle._counter;
  _id           = particle._id;
  _point_type   = particle._point_type;
  _parent_id    = particle._parent_id;
  _processor_id = particle._processor_id;
  _elem_id      = particle._elem_id;
  _neighbor_list= particle._neighbor_list;
  _force        = particle._force;
  _orientation  = particle._orientation;
}


// ======================================================================
PointParticle::~PointParticle()
{
  // do nothing
}

  
  
// ======================================================================
void PointParticle::set_particle_force(const std::vector<Real>& pforce)
{
  _force = pforce;
}
  
  

  
// ======================================================================
void PointParticle::add_particle_force(const std::vector<Real>& pforce)
{
  for(std::size_t i=0; i<pforce.size(); ++i)
  {
    _force[i] += pforce[i];
  }
}


  
// ======================================================================
void PointParticle::zero_particle_force()
{
  for(std::size_t i=0; i<3; ++i)
  {
    _force[i] = 0.0;
  }
}
  
  

// ======================================================================
void PointParticle::reinit_particle()
{
  // reset int
  _processor_id = -1;
  _elem_id      = -1;

  // reinit vectors
  this->zero_particle_force();

  _neighbor_list.clear();

  // reinit orientation vector
//  for(std::size_t i=0; i<_orientation.size(); ++i){
//    _orientation[i] = 0.0;
//  }
  
}
  


// ======================================================================
void PointParticle::set_orientation(const std::vector<Real>& rot_vec)
{
  _orientation = rot_vec;
}



// ======================================================================
void PointParticle::print_info(const bool & print_neighbor_list) const
{
  /* --------------------------------------------------------------------------
   * Scheme 1: using printf. Then every process will print out info on its own.
   * --------------------------------------------------------------------------*/
  printf("point particle[%d]: \n", _id);
  printf("      center = (%f, %f, %f)\n", _center(0), _center(1), _center(2));
  printf("      PBC counter = (%d, %d, %d)\n", _counter[0],_counter[1],_counter[2]);
  printf("      force  = (%f, %f, %f)\n", _force[0],_force[1],_force[2]);
  
  // output elem id and process id
  printf("      parent_id   = %d\n",    _parent_id);
  printf("      point_type  = %d\n",   _point_type);
  printf("      elem_id     = %d  \n",    _elem_id);
  printf("      proc_id     = %d\n", _processor_id);
  
  // output orientation vector
  printf("      orientation = (" );
  for(std::size_t i=0; i<_orientation.size(); ++i)
  {
    printf("%f",_orientation[i]);
    if(i != _orientation.size()-1 )printf(", ");
  }
  printf(")\n");
  
  // output the neighbor list
  if (print_neighbor_list)
  {
    if( _neighbor_list.size()>0 )
    {
      printf("      its neighbor list includes(%lu): \n",_neighbor_list.size());
      for (std::size_t i=0; i<_neighbor_list.size(); ++i )
        printf("      ---point %lu,   distance = %f\n",
               _neighbor_list[i].first, _neighbor_list[i].second);
    }
    else
      printf("      there are no neighbors around this point particle!\n");
  }
  printf("\n");
  
} // the end of print_info()

  
  
}   // end of namespace
