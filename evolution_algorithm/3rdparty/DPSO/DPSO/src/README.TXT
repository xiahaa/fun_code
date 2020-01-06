@ DPSO is based on the PSO version by Yuhui SHI. 
@ Please see the OLD_README before reading this file!!!

/====================  Other Information =============================/
  Dissipative PSO V1.1  
  Xiao-Feng Xie, Dec 15, 2001 (last modified)
  URL: http://www.adaptivebox.net/research/

  Difference: add two chaos parameters: Cv, Cl
  @ Ref: X. F. Xie, W. J. Zhang, Z. L. Yang. A dissipative particle swarm optimization.
    Congress on Evolutionary Computation (CEC), Hawaii, USA, 2002: 1456-1461 

  @@Notes: IF Cv=Cl=0, then DPSO==PSO

/========================= Usage ============================
To make all the programs you should type 
->	make all

To run it, under command prompt, type in
->	dpso  *.run 

@ The *.run file specifies all information required to run 
the program. There are (you may see the C Program for dedails): 

/*------Algorithm setting parameters: Total evaluation time = m*(T+1)------------*/
NUMBER_OF_AGENTS: m
MAXITER: T
weight_up
weight_low
Cv
Cl
MAXV

/*------Output setting parameters------------*/
run_no
Output step of generation
results_file_prefix

/*------problem setting parameters------------*/
fun_type: 0-Shaffer f6, 1-Sphere, 2-Rosenbrock, 3-generalized Rastrigrin, 4-generalized Griewank
DIMENSION
MAXX

/========================= Examples ============================
* examples in the directory "testfile", All are under m=20, and Dim=10
SF1: SPSO, Inertia weight(0.9 to 0.4)
Griewank_LPSO_SF1.RUN
Rastrigrin_LPSO_SF1.RUN

SF3: DPSO, Inertia weight(0.4), cv=0, cl=0.001
Griewank_DPSO_DP3.RUN
Rastrigrin_DPSO_DP3.RUN

/========================= Output Information ============================
The output file name will be given in runtime.

Each row represents the information of a run. hence the total row number is equal to the parameter "run_no" if it was not interrupted.
Each column represents the gbest value at a certain generation number, which the interval is equatio to the parameter "Output step of generation".
