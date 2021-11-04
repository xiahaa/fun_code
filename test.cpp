#include "mathbasic.h"
#include <stdio.h>


#define M_PI 3.1415926535897932384626433 
#define DEG2RAD M_PI/180.0
#define RAD2DEG 180.0/M_PI



int main(int argc, char **argv)

{
	float eulers[3] = {0};
	printf("Enter Pitch Yaw Roll angles\n");
	scanf("%f %f %f",&eulers[0],eulers+1,eulers+2);
	printf("Euler Angles: Roll %f ; Pitch %f ; Yaw %f \n",eulers[2],eulers[1],eulers[0]);
	eulers[0] = eulers[0] * DEG2RAD;
	eulers[1] = eulers[1] * DEG2RAD;
	eulers[2] = eulers[2] * DEG2RAD;

	float r[9];	
	e2r(eulers,r);
	printf("e2r: %f ; %f ; %f ;\n %f ; %f ; %f ;\n %f ; %f ; %f \n",r[0],r[1],r[2],
																	r[3],r[4],r[5],
																	r[6],r[7],r[8]);
	float e[3];
	r2e(r,e);
	printf("q2e: R %f ; P %f; Y %f \n",e[2] * RAD2DEG,e[1]* RAD2DEG,e[0]* RAD2DEG);



	int wait;
	scanf("%d",&wait);

	return 0;

}