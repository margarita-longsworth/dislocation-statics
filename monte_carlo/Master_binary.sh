#! /bin/bash

#mpic++ -fopenmp -std=c++0x -o MasterSlaveHybrid MasterSlaveHybrid.c++ Particle.o Cell.o Domain.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o -L/home/users/ganeshfx/bin -limd
#mpic++ -std=c++0x -o MasterSlave MasterSlave.c++ Particle.o Cell.o Domain.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o -L/home/users/ganeshfx/bin -limd

#mpic++ -std=c++0x -o MasterSlavePrecheck MasterSlavePrecheck.c++ Particle.o Cell.o Domain.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o -L/home/users/ganeshfx/bin -limd

mpic++ -std=c++0x -o MasterSlavePostTest MasterSlavePostTest.c++ Particle.o Cell.o Domain.o localMD.o acceptanceCheck.o TrialMove.o UtilityMethods.o -L/home/users/ganeshfx/bin -limd
