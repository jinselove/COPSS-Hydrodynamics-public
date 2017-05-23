#include "copss_point_particle_system.h"

// this is just for test
int main (int argc, char**argv){
	CopssPointParticleSystem system(argc, argv);

	time_t rawtime;
  	struct tm * timeinfo;
  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );	
	system.start_time(timeinfo);
	system.init_system("polymer_control.in");
	system.end_time(timeinfo);
	return 0;
}