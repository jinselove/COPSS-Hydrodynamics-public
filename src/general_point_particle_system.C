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

// void GeneralPointParticleSystem::read_data(std::string control_file){
//     const GetPot tmp(control_file);
//     input_file = tmp;
// 	this -> read_test_name();
// 	this -> read_physical_parameter();
// 	this -> read_particle_parameter();
// }


void GeneralPointParticleSystem::read_particle_parameter(){

	const string particle_type = input_file("particle_type","other");
	if(particle_type != "point_particle" and particle_type != "finite_size_particle"){
		error_msg = "	Invalid particle_type !!!";
		libmesh_error();    
	}       
	point_particle_model	= input_file("point_particle_model", "other");
	if(point_particle_model == "bead"){
		Nb = input_file("Nb",21); // total # of point particles
		Ns = Nb - 1;//we can look at these beads as a single chain without spring force. Just for development convinence
	}
	else if(point_particle_model == "polymer_chain"){
	Nb = input_file("Nb", 21);// total # of beads
		Ns = input_file("Ns", 20);// # of springs per Chain
		nChains = Nb / (Ns+1);
		bk = input_file("bk", 1E-6);//Kuhn length(um)
		Nks = input_file("Nks", 1E-6);// # of Kuhn length per spring
		Ss2 = Nks*bk*bk/6.; // (um^2)
		q0 = Nks * bk;// maximum spring length (um)
		chain_length = Ns * q0; // contour length of the spring (um)
		Dc = Db / Real(Nb); // Diffusivity of the chain (um^2/s)
	}
	else{
	error_msg = "	Invalid point_particle_model !!!";
	PMToolBox::output_message(error_msg, comm_in);
	libmesh_error();
	}
	//make sure we have the right combination of Nb and Ns
	if((Nb % (Ns + 1)) != 0){
		error_msg = "	Incorrect combination of Nb (number of beads) and Ns (number of springs per chain)";
		PMToolBox::output_message(error_msg, comm_in);
		libmesh_error();     
	}

	if(comm_in.rank() == 0){
		printf("##########################################################\n"
		       "#                  Particle Parameters                    \n"
		       "##########################################################\n\n"
		       "   Particle type             : %s\n"
		  	   "   Point Particle model      : %s\n"
		 	   "   number of beads       Nb  = %d\n",
		 	   particle_type.c_str(), point_particle_model.c_str(), Nb);

		  // for particular models
		if(point_particle_model == "polymer_chain"){
			printf( "   number of springs per Chain       Ns  = %d\n"
					"   number of Chains              nChains = %d\n"
					"   Kuhn length                       bk  = %.6e (um)\n"
					"   # of Kuhn segment per spring      Nks = %.6e\n"
					"   second moment of polymer chain    Ss2 = %.6e (um^2)\n"
					"   maximum spring length             q0  = %.6e (um)\n"
					"   chain length of polymer           Lc  = %.6e (um)\n"
					"   chain diffusivity                 Dc  = %.6e (um^2/s)\n",
		  			Ns, nChains, bk, Nks, Ss2, q0, chain_length, Dc);
		}

		printf("------------> The non-dimensional variables:\n");


		printf( "   non-dimensional bead radius      a0     = %.1e\n"
			  	"   non-dimensional ksi = sqrt(PI)/(3a0)    = %.6e\n",
			  	1.0, std::sqrt(PI)/(3.) );
		if(point_particle_model == "polymer_chain"){
			printf( "   non-dimensional Kuhn length    bk/a     = %.6e\n"
			  		"   non-dimensional spring length  q0/a     = %.6e\n"
			  		"   non-dimensional contour length Lc/a     = %.6e\n"
			  		"   non-dimensional Ss/a = sqrt(Ss2/a^2)    = %.6e\n"
			  		"   non-dimensional ksi = sqrt(PI)/(3a0)    = %.6e\n",
			  		bk/Rb, q0/Rb, chain_length/Rb, std::sqrt(Ss2/Rb/Rb), std::sqrt(PI)/(3.) );
		}
	} // end if (comm_in.rank() == 0)
}

} // end of namespace