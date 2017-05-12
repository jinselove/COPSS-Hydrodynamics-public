#include "general_point_particle_system.h"

using std::cout;
using std::endl;
using std::string;

namespace libMesh{


GeneralPointParticleSystem::GeneralPointParticleSystem(int argc, char** argv)
	:Copss(argc, argv)
{
	// do nothing
}


void GeneralPointParticleSystem::build_system(){
	this -> read_data("polymer_control.in");
}


}