/*
  This program optimizes benchmark functions using particle swarm algorithm
  Asynchronous version
  Yuhui Shi, May 15, 1998
  @ See OLD_README file for more details
  
  /-----------------------------------------------------------------------------------
  Dissipative PSO V1.1
  Xiao-Feng Xie, Dec 15, 2001 (last modified)
  URL: http://www.adaptivebox.net/research/
  @ See README.txt file for more details

  Difference: add two chaos parameters: Cv, Cl
  @ Ref: X. F. Xie, W. J. Zhang, Z. L. Yang. A dissipative particle swarm optimization.
    Congress on Evolutionary Computation (CEC), Hawaii, USA, 2002: 1456-1461 
  /-----------------------------------------------------------------------------------
*/
//normal file
#include "headfile.h"
#include "global.h"
#include "mem_loc.h"
#include "myfun.h"
#include "methods.h"

float calcFunc(int fun_type, FMATRIX outxx, int a, int DIMENSION);
float getRandom();

/* ******************** main() ******************** */
int main (int argc, char *argv[])
{
  /*------Algorithm setting parameters: Total evaluation time = m*(T+1)------------*/
  int  NUMBER_OF_AGENTS = 20;           //m: the number of particles
  int  MAXITER = 2000;                  //T: the maximum iteration cycles

  float weight_up=0.9, weight_low=0.4;  //inertia weight, for linear decreasing PSO (LPSO)
  float c1=2, c2=2;                     //learing factors
  float Cv = 0.0, Cl = 0.001;           //chaos factors: by Xiao-Feng Xie
  float MAXV = 0.5;			//For controllong the maximum velocity as MAXV*(dimension length)

  /*------Output setting parameters------------*/
  int run_no;                           //number of runs
  int OutputStep = 1;                   //For output interval
  char resfilePrefix[120];              //The prefix of output file

  /*------problem setting parameters------------*/
  int fun_type = 3;			// 0:Schaffer f6  1:sphere  2:Rosenbrock  3:generalized Rastrigrin  4:generalized Griewank
  int DIMENSION = 30;			// the dimension of problem
  float MAXX = 600;			// problem range: [-MAXX, MAXX] for all dimensions

  /*------temp values ------------*/
  FILE *fp, *frun;
  char runfile[120], resfile[120];      //run file and log file
  char temp[120];

  int a,b;
  int i,j;

  float minval=0.0;
  float realWeight;
  int iter;
  int gbest;
  int firsttime;
  int finish;


  time_t tt;
  struct tm *newtime;
  int curSec = 0;

  /* *********************************************************
		Open problem file
	********************************************************* */
	if (argc<2)
	{
		printf("Need to specify runfile. For Example: dpso *.run");
		exit(1);
	}
	strcpy(runfile,argv[1]);

	if ((frun=fopen(runfile,"r"))==NULL)
	{
		printf("Cant read file");
		exit(1);
	}

	/*Get parameters from file */
	fscanf(frun, "%d %d %f %f %f %f %f %d %d %s %d %d %f",
	  &NUMBER_OF_AGENTS, &MAXITER, &weight_up, &weight_low, &Cv, &Cl, &MAXV,
	  &run_no, &OutputStep, resfilePrefix,
	  &fun_type, &DIMENSION, &MAXX);
	fclose(frun);
	
	/*allocate memory*/
	FVectorAllocate(&trail,	 DIMENSION);
	FVectorAllocate(&pbest,  NUMBER_OF_AGENTS);
	FMatrixAllocate(&maxx,   DIMENSION, 2);
	FMatrixAllocate(&xx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&vx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&tx,     DIMENSION,NUMBER_OF_AGENTS);
	FMatrixAllocate(&pbestx, DIMENSION,NUMBER_OF_AGENTS);

	/****************************initial parameter range********************************/
	for (a=0;a<DIMENSION;a++)
	{
		maxx[a][0]=-MAXX;			  /* range of xx[]  */
		maxx[a][1]= MAXX;             /* range of xx[]  */
	}


	time(&tt);
	printf("begin time: %s\n",ctime(&tt));

	/*** .log file for output best agent index vs iteration ***/
	sprintf(temp, "(N%d, V%.2f, W%.2f-%.2f, G%d, Cv%.4f, Cl%.4f, D%d)", NUMBER_OF_AGENTS, MAXV, weight_up, weight_low, MAXITER, Cv, Cl, DIMENSION);
	strcpy(resfile,resfilePrefix);
	strcat(resfile,temp);
	strcat(resfile,".log");
	printf("Out: %s\n", resfile);
	if ((fp=fopen(resfile,"w"))==NULL)
	{
		printf("Cant write file");
		exit(1);
	}

	/* ----- loop for runs -----*/
	for (i=0;i<run_no;i++)
	{
		firsttime=1;           //first iteration of this run
		iter=0;
		gbest=0;               //initialy assume the first particle as the gbest

		time(&tt);
		newtime = localtime( &tt ); /* Convert to local time. */
		curSec = newtime->tm_sec;
		for (j=0; j< curSec; j++) {
			rand();
		}
		
  		/* **********************************************
			This loop initializes the individual agents  for each run
		********************************************** */
		for (a=0;a<NUMBER_OF_AGENTS;a++) {
			for (b=0;b<DIMENSION;b++)	{
				xx[b][a] = (float) ((maxx[b][1] - maxx[b][0])*getRandom() + maxx[b][0]);
				pbestx[b][a]=xx[b][a];
				vx[b][a] = (float) (MAXV*(maxx[b][1]-maxx[b][0])*getRandom());
				if (getRandom()> 0.5) vx[b][a]=-vx[b][a];
			}
		}
		/* *******************************************************
			Main Work Loop for each run here
		******************************************************** */

		finish=0;
		do {
			/*------ more randomize by Xiao-feng Xie start ------*/
			time(&tt);
			newtime = localtime( &tt ); /* Convert to local time. */
			curSec = newtime->tm_sec;
			for (j=0; j< curSec%7; j++) {
				rand();
			}
			/*------ more randomize by Xiao-feng Xie end   ------*/

			iter++;

			//update inertia weight
			realWeight = (weight_up-weight_low) * (MAXITER - iter) /MAXITER +weight_low;

			for (a=0;a<NUMBER_OF_AGENTS;a++)
			{
				/* *********************************************
				eval(a) error is returned by function routines 
				********************************************* */
				minval = calcFunc(fun_type, xx, a, DIMENSION);
				if (firsttime==1) pbest[a]=minval;

				if (minval < pbest[a]) {
					pbest[a]=minval;
					for (b=0;b<DIMENSION;b++) pbestx[b][a]=xx[b][a];
					if (pbest[a] < pbest[gbest]) {
						gbest=a;
					}
				}

			}

			for (a=0;a<NUMBER_OF_AGENTS;a++) {

				/* asynchronous version */
				for (b=0;b<DIMENSION;b++)	{
					vx[b][a] = (float) (realWeight*vx[b][a] + c1*getRandom()*(pbestx[b][a]-xx[b][a]) +
					c2*getRandom()*(pbestx[b][gbest]-xx[b][a]));
					if (vx[b][a]>MAXV*maxx[b][1]) {
						vx[b][a]=MAXV*maxx[b][1];
					} else if (vx[b][a]<MAXV*maxx[b][0]) {
						vx[b][a]=MAXV*maxx[b][0];
					}
					tx[b][a]=xx[b][a]+vx[b][a];
				}
			}

			for (a=0;a<NUMBER_OF_AGENTS;a++) {
				for (b=0;b<DIMENSION;b++)	{
					xx[b][a] = tx[b][a];
                    /*------ chaos factors Cv, Cl for absentminded birds: by Xiao-Feng Xie start ------*/
                    if (getRandom()< Cv) {		
						vx[b][a] = (float) (MAXV*(maxx[b][1]-maxx[b][0])*getRandom());
						if (getRandom() > 0.5) vx[b][a]=-vx[b][a];
					}
					if (getRandom()< Cl) {		
						xx[b][a] = (float) ((maxx[b][1]-maxx[b][0])*getRandom()+maxx[b][0]);
					}
                    /*------ chaos factors Cv, Cl for absentminded birds: by Xiao-Feng Xie end  ------*/
				}
			}

			/* ******************************************************
				Terminate on criterion
			********************************************************* */
			//output best index vs iteration
			if (iter%OutputStep==0) {
				if (fun_type==0) {
					fprintf(fp,"%e\t",1.0-pbest[gbest]);
				}
				else {
					fprintf(fp,"%e\t",pbest[gbest]);
				}
			}

			if (iter >= MAXITER)
			{
				fprintf(fp,"\n");
				printf("%d run: gbest=%.4f\n",(i+1), pbest[gbest]);
				finish=1;
			}
			firsttime=0;

		}     /* **************** End of do-loop *************** */
		while (! finish);
	}
	time(&tt);
	fclose(fp);
	printf("end time: %s\n",ctime(&tt));
	return 0;
}

float getRandom() {
  return rand()/(float)RAND_MAX;
}

float calcFunc(int fun_type, FMATRIX outxx, int a, int DIMENSION) {
	float minval;
	switch (fun_type)
	{
		case 0:
			minval=f6(outxx, a);
			break;
		case 1:
			minval=sphere(outxx, a,DIMENSION);
			break;
		case 2:
			minval=rosenbrock(outxx, a,DIMENSION);
			break;
		case 3:
			minval=rastrigrin(outxx, a,DIMENSION);
			break;
		case 4:
			minval=griewank(outxx, a,DIMENSION);
			break;
		default:
			printf("\n Not a valid function type\n");
			exit(1);
	}
	return minval;
}

