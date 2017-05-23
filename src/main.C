#include "general_point_particle_system.h"

// this is just for test
int main (int argc, char**argv){
	GeneralPointParticleSystem system(argc, argv);

	time_t rawtime;
  	struct tm * timeinfo;
  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );	
	system.start_time(timeinfo);
	system.build_system();
	system.end_time(timeinfo);
	return 0;
}