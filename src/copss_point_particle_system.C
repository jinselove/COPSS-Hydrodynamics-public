#include "copss_point_particle_system.h"

using std::cout;
using std::endl;
using std::string;

namespace libMesh{



//==========================================================================
CopssPointParticleSystem::CopssPointParticleSystem(CopssInit& init)
:Copss(init)
{

}

CopssPointParticleSystem::~CopssPointParticleSystem(){
	delete polymer_chain;
	delete point_mesh;
	polymer_chain = NULL;
	point_mesh = NULL;
}


//==========================================================================
void CopssPointParticleSystem::read_particle_info(){
	if (particle_type != "point_particle"){
		error_msg = "invalid particle type ("+particle_type+") defined\n";
		PMToolBox::output_message(error_msg,comm_in);
		libmesh_error();
	}
	point_particle_model	= input_file("point_particle_model", "other");
	if(point_particle_model == "bead"){
		Nb = input_file("Nb",21); // total # of point particles
		Ns = Nb - 1;//we can look at these beads as a single chain without spring force. Just for development convinence
		nBonds = Ns;
	}
	else if(point_particle_model == "polymer_chain"){
		Nb = input_file("Nb", 21);// total # of beads
		Ns = input_file("Ns", 20);// # of springs per Chain
		nChains = Nb / (Ns+1);
		nBonds = nChains * Ns;
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
}// end read_particle_parameter()


//==========================================================================
void CopssPointParticleSystem::create_object(){
  const unsigned int chain_id = 0;
  polymer_chain = new PolymerChain (chain_id, *pm_periodic_boundary);
  std::ostringstream pfilename;
  if(restart)
  {
	pfilename << point_particle_model<<"_data_restart_"<< restart_step << ".vtk";
	output_msg = "-------------> read "+point_particle_model+" data from "+pfilename.str()+ " in restart mode\n";
	PMToolBox::output_message(output_msg, comm_in);	
    polymer_chain->read_data_vtk(pfilename.str());
  } 
  else
  {
	pfilename << "point_particle_data.in";
	polymer_chain->read_data_pizza(pfilename.str(), Nb, nBonds, comm_in.rank());
	output_msg = "--------------> polymer_chain object is built for copss_point_particle_system using data from "+pfilename.str();
	PMToolBox::output_message(output_msg, comm_in);
	//comm_in.barrier();
  }
  pfilename.str(""); pfilename.clear();
}//end function

//=====================================================================
void CopssPointParticleSystem::create_object_mesh(){
  // prepare domain and objects
  this -> create_domain_mesh();
  this -> create_periodic_boundary();
  this -> create_object();

  // create object mesh
  const Real search_radius_p = 4.0/alpha;
  const Real search_radius_e = 0.5*max_mesh_size + search_radius_p;

  point_mesh = new PointMesh<3> (*mesh, *polymer_chain, search_radius_p, search_radius_e);

  point_mesh->add_periodic_boundary(*pm_periodic_boundary);

  point_mesh->reinit();

  if(comm_in.rank() == 0){
  	printf("-------------> Reinit point mesh object, finished! \n"
  		   "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
		   "### The point-mesh info:\n"
		   "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
		   "Total number of point particles: %d\n "
		   "search_radius_p = %.4e, , search_radius_e = %.4e\n\n",
		   point_mesh->num_particles(), search_radius_p, search_radius_e);
	point_mesh->print_point_info();
   }
} // end function



} // end of namespace
