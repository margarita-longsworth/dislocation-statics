#! /bin/bash

rm *.a *.o


mpic++ -c -std=c++0x -o Particle.o Particle.c++
mpic++ -c -std=c++0x -o Cell.o Cell.c++ 
mpic++ -c -std=c++0x -o Domain.o Domain.c++
#mpic++ -c -std=c++0x -o constructSampleZone.o constructSampleZone.c++
mpic++ -c -std=c++0x -o localMD.o localMD.c++
mpic++ -c -std=c++0x -o UtilityMethods.o UtilityMethods.c++
mpic++ -c -std=c++0x -o acceptanceCheck.o acceptanceCheck.c++
mpic++ -c -std=c++0x -o TrialMove.o TrialMove.c++

#mpic++ -c -std=c++0x -o MonteCarloReader.o MonteCarloReader.c++
#mpic++ -c -std=c++0x -o MonteCarloAPI.o MonteCarloAPI.c++


# dynamic so 
#mpic++ -fPIC -shared -o libMontecarloAPI.so *.o 


# static archive
#ar -rcs libMontecarloAPI.a Particle.o Cell.o Domain.o constructSampleZone.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o MonteCarloAPI.o

ar -rcs libMontecarloAPI.a Particle.o Cell.o Domain.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o MonteCarloReader.o

ar t libMontecarloAPI.a

echo *.o
