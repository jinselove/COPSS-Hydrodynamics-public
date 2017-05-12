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

#include "copss.h"

using std::cout;
using std::endl;
using std::string;

namespace libMesh{

class GeneralPointParticleSystem : public Copss
{
public:
	
	GeneralPointParticleSystem (int argc, char** argv);

	void build_system();

protected:
	std::string _point_particle_model;
	unsigned int _Nb; // total # of beads
	unsigned int _Ns = _Nb - 1; // total # of springs per Chain --> _Ns = _Nb-1
	unsigned int _nChains; // total # of chains
	Real _bk; // Kuhn length (um)
	Real _Nks; // # of Kuhn length per Chain
	Real _q0 = _Nks * _bk; //maximum spring length (um)
	Real _chain_length = _Ns * _q0; // contour length of the spring (um)
	
};


}
