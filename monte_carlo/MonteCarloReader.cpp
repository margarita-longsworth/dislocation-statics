/*
 * MonteCarloReader.c++
 *
 */
#include <iostream>
#include "MonteCarloReader.h"

using namespace std;

void readMCParamFile(string paramFile){

    ifstream paramin(paramFile,std::ios_base::in);
    string keyword;

    while (!paramin.eof()) {

        if (paramin.get() != '#'){

              paramin >> keyword;
              if(keyword == "cpuDimension"){
                      paramin >> mccpuDimension.x;
                      paramin >> mccpuDimension.y;
                      paramin >> mccpuDimension.z;
              }
              else if(keyword == "MCSteps"){
                     paramin >>  mcnSteps;
              }

              else if(keyword == "totalTypes"){
                     paramin >>  mcTotalTypes;
              }

              else if(keyword == "realTypes"){
                     paramin >>  mcRealTypes;
              }

              else if(keyword == "SimboxX"){
                     paramin >>  mcSimboxX.x;
                     paramin >>  mcSimboxX.y;
                     paramin >>  mcSimboxX.z;
              }

              else if(keyword == "SimboxY"){
                     paramin >>  mcSimboxY.x;
                     paramin >>  mcSimboxY.y;
                     paramin >>  mcSimboxY.z;
              }

              else if(keyword == "SimboxZ"){
                     paramin >>  mcSimboxZ.x;
                     paramin >>  mcSimboxZ.y;
                     paramin >>  mcSimboxZ.z;
              }

              else if(keyword == "mczoneSize"){
                     paramin >>  mcZoneSize.x;
                     paramin >>  mcZoneSize.y;
                     paramin >>  mcZoneSize.z;
              }

              else if(keyword == "pbc"){
                      paramin >> mcPBC.x;
                      paramin >> mcPBC.y;
                      paramin >> mcPBC.z;
              }

              else if(keyword == "temperature"){
                      paramin >> mcTemperature;
              }
              else if(keyword == "sphereRadius"){
                      paramin >> mcSphereRadius;
              }

              else if(keyword == "sphereWall"){
                      paramin >> mcSphereWallThickness;
              }

              else if(keyword == "cohesiveEnergy"){
                      paramin >> mcCohesive;
              }

              else if(keyword == "flushInterval"){
                      paramin >> mcFlushInterval;
              }

              else if(keyword == "mdInput"){
                      paramin >> mcInputFile;
              }

              else if(keyword == "sphereParam"){
                      paramin >> mcSphereFile;
              }


           }

    }

    cout << " Parameter File reading Successfull!!!" << endl;

    cout << " input total Steps   : "  << mcnSteps << endl;
    cout << " input CPU Dimension : [ "<< mccpuDimension.x <<" "<<mccpuDimension.y<<" "<<mccpuDimension.z<<" ]"<< endl;
    cout << " Simbox X            : [ "<< mcSimboxX.x <<" "<< mcSimboxX.y<<" "<<mcSimboxX.z<<" ]"<< endl;
    cout << " Simbox Y            : [ "<< mcSimboxY.x <<" "<< mcSimboxY.y<<" "<<mcSimboxY.z<<" ]"<< endl;
    cout << " Simbox Z            : [ "<< mcSimboxZ.x <<" "<< mcSimboxZ.y<<" "<<mcSimboxZ.z<<" ]"<< endl;
    cout << " Zone Size           : [ "<< mcZoneSize.x <<" "<< mcZoneSize.y<<" "<<mcZoneSize.z<<" ]"<< endl;
    cout << " Periodic Boundary   : [ "<< mcPBC.x <<" "<< mcPBC.y<<" "<<mcPBC.z<<" ]"<< endl;
    cout << " input temperature   : "  << mcTemperature << endl;

    cout << " Sphere radius       : "  << mcSphereRadius << endl;
    cout << " Sphere Wall         : "  << mcSphereWallThickness << endl;
    cout << " Cohesive Energy     : "  << mcCohesive << endl;

    cout << " Flush Interval      : "  << mcFlushInterval << endl;
    cout << " mdInput file        : "  << mcInputFile << endl;
    cout << " sphere Param file   : "  << mcSphereFile << endl;

}

void readMDConfiguration(string fileName){

	// configuration file reader
	const int LIMIT=30000;
    cout<< " Input Configuration File : " << fileName << endl;

    ifstream fin(fileName,std::ios_base::in);

    // local value holders from file
    long number; int type;
    double mass,epot,eamRho;
    Vector3d position,velocity;

    char headerline[LIMIT];
    cout << " File header check " << endl;

    // crunch header part
    while (fin.get() == '#'){
          fin.getline(headerline,LIMIT,'\n');
          cout << headerline << endl;
          continue;
    }

    // read and fill MC vector container

    while((!fin.eof())){ // end of file check
    	   // reading            // pushing
           fin>>number;        mcNumbers.push_back(number);
           fin>>type;          mcTypes.push_back(type);
           fin>>mass;          mcMasses.push_back(mass);
           fin>>position.x;    mcPositions.push_back(position.x);
           fin>>position.y;    mcPositions.push_back(position.y);
           fin>>position.z;    mcPositions.push_back(position.z);
           fin>>velocity.x;    fin>>velocity.y; fin>>velocity.z;
           fin>>epot;          mcPotentialEnergies.push_back(epot);
           fin>>eamRho;
           fin.get(); //Linebreak
           fin.peek(); //Check for eof
    }

    fin.clear();      // resetting bit states
    fin.close();

    cout << "  ======= Inside File reader-- container Check  =========== " << endl;
    cout << " mcNumbers           size : " << mcNumbers.size() << endl;
    cout << " mcTypes             size : " << mcTypes.size() << endl;
    cout << " mcMasses            size : " << mcMasses.size() << endl;
    cout << " mcPositions         size : " << mcPositions.size()/3 << endl;
    cout << " mcPotentialEnergies size : " << mcPotentialEnergies.size() << endl;
}


