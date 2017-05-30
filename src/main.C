#include "copss_point_particle_system.h"

// this is just for test
int main (int argc, char **argv){

  CopssInit init(argc, argv);


  CopssPointParticleSystem system(init);

  
  // Print out the current date and time 
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
//	system.start_time(timeinfo);

	system.init_system("polymer_control.in");
  
	return 0;
}