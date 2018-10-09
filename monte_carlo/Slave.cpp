/*
 * Slave.c++
 *
 * Slave Class contains an important method named runAsSlave which works in coordination with the Master Class with MPI communications.
 * The prime role of the runAsSlave method is
 *  1. to receive Job from Master
 *  2. to execute a local MD
 *  3. to perform acceptance Check
 *  4. export back configuration to Master (if accepted)
 */

#include "Slave.h"

using namespace std;

Slave::Slave(int MasterID,string fileName,MPI_Comm subComm){
	 MasterID_    = MasterID;
	 fileName_    = fileName;
	 subComm_     = subComm;
}

Slave::~Slave(){

}

void Slave::runAsSlave(void){

	const int TAG_RUNJOB=0, TAG_TERMINATE=1, TAG_FREE=2;
	long available = 1;

	// Post process availability
	MPI_Send(&available, 1, MPI_LONG, MasterID_, TAG_FREE, MPI_COMM_WORLD);

    int slaveRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&slaveRank);

    interfaceIMD        mdObject;
    Particle            refParticle;
	vector<Particle>    inputSphere;
	vector<Particle>    updatedSphere;
	vector<Particle> testConfig; // debugging

	// getting char* from string
	char *paramFile = new char[fileName_.length() + 1];

    std::strcpy(paramFile, fileName_.c_str());
	mdObject.loadInputFile(paramFile);

	bool slaveFree=true;
	double workerResult[5];

	//  Sphere Core radii + SphereWall Thickness + SphereWall Thickness;
    string outFile;
    int convergenceFlag=0;

	while(slaveFree){

		long buffInd;
		long particleCounter=0;
		long bufferInfo[5];

		// Receive buffer size
	    MPI_Recv(bufferInfo, 5, MPI_LONG, MasterID_, MPI_ANY_TAG,
		       MPI_COMM_WORLD, &status_);

	    // buffer Size
	    long bufferLength = bufferInfo[0];
        long stepNo       = bufferInfo[1];
	    long versionNo    = bufferInfo[2];
	    long rtoVID       = bufferInfo[3];  // ID of particle which is to be transmuted from real to virtual
        long vtoRID       = bufferInfo[4];  // ID of particle which is to be transmuted from virtual to real

        bool is_finTemp;

        double eDiff;

        //energyDiff = &eDiff;

//        cout<<"///////////////////////////////////"<<endl;
//        cout<<" Inside Slave check for values"<<endl;
//        cout<<" rtoVID : "<<rtoVID<<endl;
//        cout<<" vtoRID : "<<vtoRID<<endl;
//        cout<<"///////////////////////////////////"<<endl;

        if((rtoVID==-1) && (vtoRID==-1)){ // T=0K
        	is_finTemp = false;
        }
        else{
        	is_finTemp = true;
        }

	    switch (status_.MPI_TAG)
		{
	        case TAG_RUNJOB:
		       {
			     // Allocate buffers
			     double* localBuffer;
			     localBuffer = new double [bufferLength];

		    	 // Receive Buffers
		    	 MPI_Recv(localBuffer, bufferLength, MPI_DOUBLE, MasterID_, TAG_RUNJOB,
		    		       MPI_COMM_WORLD, &status_);

		    	 particleCounter = bufferLength/7; // each particle has 7 attributes
		    	 buffInd         = 0;

		    	 // Fill particle container from buffer
		    	 for(auto i=0; i<particleCounter; i++){

		    		 Particle  particleObject;
		    		 particleObject.setmcID(static_cast<long> (localBuffer[buffInd++]));
                     particleObject.setType(static_cast<int> (localBuffer[buffInd++]));
                     particleObject.setMass(localBuffer[buffInd++]);

                     double posX = (localBuffer[buffInd++]);
                     double posY = (localBuffer[buffInd++]);
                     double posZ = (localBuffer[buffInd++]);

                     particleObject.setPosition(posX,posY,posZ);

                     // for debugging
                     if(particleObject.getmcID()==rtoVID){
                    	 Vector3d var_pos = particleObject.getPosition();
                    	 cout<<"Chosen carbon coordinate :"<<var_pos.x<<" "<<var_pos.y<<" "<<var_pos.z<<endl;
                     }
                     if(particleObject.getmcID()==vtoRID){
                    	 Vector3d var_pos = particleObject.getPosition();
                    	 cout<<"Chosen vatom coordinate :"<<var_pos.x<<" "<<var_pos.y<<" "<<var_pos.z<<endl;
                     }

                     particleObject.setEpot(localBuffer[buffInd++]);
                     inputSphere.push_back(particleObject);
		    	 }

		    	 // invoke Local MD call
		    	 // for T=0K;  convergence flag = '1' if simulation converged otherwise '0'
		    	 // for T>0K;  always convergence flag = '1'
                 //cout << "Slave Processor "<<slaveRank <<" is running step No  "<<stepNo <<" with "<<particleCounter
                 //		 <<" particles "<< endl;

		    	 // for debugging
            	 string infileName  = "in_config_"+to_string(stepNo)+"_"+to_string(versionNo)+".chkpt";
            	 //string outfileName = "out_config_"+to_string(stepNo)+to_string(versionNo)+".chkpt";

		    	 if(is_finTemp){
//                	 cout<<" Slave running finiteTemp MD"<<endl;
                	 //writeSphereFile(infileName,inputSphere); // write the local domain

                	 // ---- For debugging
                	 //string check_file = "in_config_369_0.chkpt";
                	 //string check_file = "in_config_639_1.chkpt";
                	 //string check_file = "in_config_818_3.chkpt";


                	 //getConfigurationFromFile(check_file,testConfig);

                	 //rtoVID=3872224; //369_0
                	 //vtoRID=498749;

//                	 rtoVID=4008331;
//                     vtoRID=1976980;

                	 //rtoVID=344935;  // 818_3
                	 //vtoRID=2288896;
                	 // ----

		    		 cout<<"IN CHECK -- Step No :"<<stepNo <<" Version no : "<< versionNo<< endl;

		    		 // convergenceFlag = mdObject.runLocalMD_finiteTemp(testConfig,updatedSphere,eDiff,rtoVID,vtoRID,subComm_); // debug hack
                	 convergenceFlag = mdObject.runLocalMD_finiteTemp(inputSphere,updatedSphere,eDiff,rtoVID,vtoRID,subComm_);
                	 cout<<"OUT CHECK -- Step No :"<<stepNo <<" Version no : "<< versionNo<< endl;

                	 std::remove(infileName.c_str()); // delete the input configuration, if MD run is successful

//                	 writeSphereFile(outfileName,updatedSphere);
//                	 cout<<" Finite temperature run successful"<<endl;
                 }
                 else{
                	 //cout <<"WARNING: Why I am coming here ? "<<endl;
                	 convergenceFlag = mdObject.runLocalMD(inputSphere,updatedSphere,subComm_);
                 }

			     // flush the intermediate Sphere Configuration
		    	 outFile  = "outputSphere_P"+to_string(slaveRank)+"_"+to_string(stepNo)+".chkpt";
		    	 //writeSphereFile(outFile,updatedSphere);

                 long exportLength  = (updatedSphere.size()*7);
                 workerResult[0] = convergenceFlag;
                 workerResult[1] = stepNo;
                 workerResult[2] = versionNo;
                 workerResult[3] = exportLength;

                 if(is_finTemp){ // finiteTemp case
                	 workerResult[4] = eDiff;
                	 //cout<<" Received energy Diff : "<<eDiff<<endl;
                 }
                 else{           // T=0K case
                	 workerResult[4] = 0;
                 }

                 MPI_Send(workerResult, 5, MPI_DOUBLE, MasterID_, TAG_RUNJOB, MPI_COMM_WORLD);
                 bool isRequired = false;
		    	 MPI_Recv(&isRequired, 1, MPI_C_BOOL, MasterID_, TAG_RUNJOB,
		    		       MPI_COMM_WORLD, &status_);

		    	 if(isRequired){
					 // collect Results
					  double* exportBuffer = new double[exportLength];
					  long expInd       = 0;
					  for(const auto& p : updatedSphere){
							 exportBuffer[expInd++] = static_cast<double> (p.getmcID());
							 exportBuffer[expInd++] = static_cast<double> (p.getType());
							 exportBuffer[expInd++] = p.getMass();
							 exportBuffer[expInd++] = p.getPosition().x;
							 exportBuffer[expInd++] = p.getPosition().y;
							 exportBuffer[expInd++] = p.getPosition().z;
							 exportBuffer[expInd++] = p.getEpot();
					  }
	                 // export buffer
	                 MPI_Send(exportBuffer,exportLength, MPI_DOUBLE, MasterID_, TAG_RUNJOB, MPI_COMM_WORLD);
	                 delete[] exportBuffer;
		    	 }

		    	 delete[] localBuffer;

		    	 inputSphere.clear();
		    	 updatedSphere.clear();
		    	 //cout << "Slave Processor "<< slaveRank<<" finished MD run for Step no "<<stepNo<<" with covergence flag "<<convergenceFlag << endl;
		         break;
		      }

		    // Message tag to terminate Slave process
		    case TAG_TERMINATE:
		      {
    		     slaveFree = false;
    		     cout << "Slave Process ID :"<<slaveRank<<" Normal termination Successful" << endl;
    		     break;
		      }
		}
	}
	delete[] paramFile;
}

void Slave::writeSphereFile(string outFile,vector<Particle>& updatedSphere){

                 //double boxSize = 500.0;
                 ofstream fout(outFile, ios_base::out);

		    	 // print IMD header -- check
		    	 fout <<"#F A 1 1 1 3 0 1"<< endl;
		    	 fout <<"#C number type mass x y z Epot"<< endl;

		    	 // edge configuration (Ruda)
//		    	 fout <<"#X "<<4.7984899999999999e+02<<" "<<0.0<<" "<<0.0<<" "<< endl;
//		    	 fout <<"#Y "<<0.0<<" "<<4.8065499999999997e+02<<" "<<0.0<<" "<< endl;
//		    	 fout <<"#Z "<<0.0<<" "<<0.0<<" "<<6.9959599999999995e+01<<" "<< endl;

		    	 // screw configuration (Ruda)
//		    	 fout <<"#X "<<4.8272199999999998e+02<<" "<<0.0<<" "<<0.0<<" "<< endl;
//		    	 fout <<"#Y "<<0.0<<" "<<4.8065499999999997e+02<<" "<<0.0<<" "<< endl;
//		    	 fout <<"#Z "<<0.0<<" "<<0.0<<" "<<6.9256500000000003e+01<<" "<< endl;

		    	 // screw configuration (Veiga)
		    	 fout <<"#X "<<4.8000000000000000e+02<<" "<<0.0<<" "<<0.0<<" "<< endl;
		    	 fout <<"#Y "<<0.0<<" "<<4.8000000000000000e+02<<" "<<0.0<<" "<< endl;
		    	 fout <<"#Z "<<0.0<<" "<<0.0<<" "<<6.9291499999999999e+01<<" "<< endl;

		    	 fout <<"#E "<< endl;

		    	 // loop over particles in the container
		    	 for(const auto& p : updatedSphere){
		    	       fout << p.getmcID()
		    	       << "  " << p.getType()
		    	       << "  " << setprecision(6) << p.getMass()
		    	       << "  " << setprecision(6) << p.getPosition().x
		    	       << "  " << setprecision(6) << p.getPosition().y
		    	       << "  " << setprecision(6) << p.getPosition().z
		    	       << "  " << setprecision(6) << p.getEpot()
		    	       << endl;
		    	 }

		    	 fout.close(); // closing outfile connection
}

void Slave::getConfigurationFromFile(string fileName,vector<Particle>& configParticles){

	// configuration file reader
	const int LIMIT=30000;
    cout<< " Input Configuration File : " << fileName << endl;

    ifstream fin(fileName,std::ios_base::in);
    if (!fin.good()){
    	cout<< " Input file does not exist!" << endl;
    	exit(-1);
    }

    // local value holders from file
    long loc_number; int loc_type;
    double loc_mass,loc_epot,loc_eamRho;
    Vector3d loc_position,loc_velocity;

    char headerline[LIMIT];
    cout << " File header check " << endl;

    // crunch header part
    while (!fin.eof()){
          fin.getline(headerline,LIMIT,'\n');
          cout << headerline << endl;
          if(headerline[1]=='E'){
            break;
          }
    }

    Particle particle;

    long num=1; // Particle ID (scope is completely local and hold no relevance to global particle ID)

    // read and fill MC vector container
    while((!fin.eof())){ // end of file check
    	   // reading            // setting
           fin>>loc_number;       particle.setmcID(loc_number);
           fin>>loc_type;         particle.setType(loc_type);
           fin>>loc_mass;         particle.setMass(loc_mass);
           fin>>loc_position.x;
//           if (mcPBC.x){
//        	   while (position.x < 0.) position.x += mcSimboxX.x;
//           	   while (position.x >= mcSimboxX.x) position.x -= mcSimboxX.x;
//           }
           fin>>loc_position.y;

//           if (mcPBC.y){
//        	   while (position.y < 0.) position.y += mcSimboxY.y;
//        	   while (position.y >= mcSimboxY.y) position.y -= mcSimboxY.y;
//           }
           fin>>loc_position.z;
//           if (mcPBC.z){
//        	   while (position.z < 0.) position.z += mcSimboxZ.z;
//        	   while (position.z >= mcSimboxZ.z) position.z -= mcSimboxZ.z;
//           }
           particle.setPosition(loc_position.x,loc_position.y,loc_position.z);
           //fin>>loc_velocity.x; fin>>loc_velocity.y; fin>>loc_velocity.z;
           fin>>loc_epot;          particle.setEpot(loc_epot);
           //fin>>loc_eamRho;

           // pushing particle into list
           configParticles.push_back(particle);

           fin.get(); //Linebreak
           fin.peek(); //Check for eof
           num++; // local Particle ID counter
    }

    fin.clear();      // resetting bit states
    fin.close();

    cout<<"Intermediate particle container Size : "<<configParticles.size()<<endl;

}
