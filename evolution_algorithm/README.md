# evolution_algorithm

try out of some so called population based optimization method:
## PSO
## JAYA
## TLBO

My personal opinion about this kind of algorithms would be:
1. the performance of PSO highly depends on the tuning: increase the phi of global best, converge fastly but high possibility of traped into a local minimum. increase the inertial, converge slowly or even no convergency...
2. JAYA sounds like a tuning free solution, however, you still need to assign the maximum iteration cnt. The worse is that no guarantee of convergency can be made, let alone the optimality...
3. TLBO seems to be very popular, almost 1300 citations since published. Tuning free except the maximum iteration. No guarantee of convergency.

May work towards some problems, even work very well for some specific problems. Mathematical explanation is unclear since it uses a lot of random selection.
