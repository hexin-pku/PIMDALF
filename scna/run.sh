#!/bin/bash
for(( i=1;i<=100;i++  ))
do
    mpirun ./dyn.x 
    mv distribution.out distribution${i}.out
    mv Broken_Traj.out Broken_Traj${i}.out
    mv MC_Acceptance.out MC_Acceptance${i}.out
done
