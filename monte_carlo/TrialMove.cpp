/*
 * TrialMove.c++
 *
 */
#include "TrialMove.h"
using namespace std;
TrialMove::TrialMove(int nTrialmoves,int inSeed){

	// Seeding number generator
    std::random_device rd;

    generatorTrial.seed(inSeed*2);
	nTrialmoves_         = nTrialmoves;
	distributionTrial    = std::uniform_int_distribution<int>(TARGET,SAMPLE);	//Either 1 or 2
}
TrialMove::~TrialMove(void){

}
int TrialMove::getNtrialmove(){
	return nTrialmoves_;
}
int TrialMove::getRandomNo(){
	int result;
	result = distributionTrial(generatorTrial);
	return result;
}

addDeleteParticle::~addDeleteParticle(void){
}

vector<Swap> addDeleteParticle::getSwapList(Domain& dom, long currentMCStep){	

	vector<Swap> swapList;

	if(getNtrialmove() == 1){
		int chosenType  = getRandomNo();
		long nTargets   = dom.getTargetListSize();
		Swap swap;
		// Random particle Selection (in Case of no available targets,revise decision)
		if((chosenType == TARGET) && (nTargets>1)){ // Target Particle deletion
			chosenType   = TARGET;
			swap.type    = SAMPLE;
		}else {      // Target Particle inclusion
			chosenType   = SAMPLE;
	         swap.type   = TARGET;
		}
		swap.ID  = dom.getRandomId(chosenType);
		swapList.push_back(swap);
	}
	else if(getNtrialmove() == 2){

		// Perform swap operation
		// swap operation

		if(dom.getTargetListSize()==0){
			cerr<<"ERROR: SWAP OPERATION NOT POSSIBLE FOR EMPTY TARGET LIST" <<endl;
            exit(1);
		}
		else if(dom.getSampleListSize()==0){
			cerr<<"ERROR: SWAP OPERATION NOT POSSIBLE FOR EMPTY SAMPLE LIST" <<endl;
            exit(1);
		}

		if ( dom.getSamplingMode()==0 ){	

			long insertID = dom.getRandomId(SAMPLE); 
			long deleteID = dom.getRandomId(TARGET);	

			Swap swap1,swap2;

			swap1.ID    = insertID;
			swap1.type  = TARGET;
			swap2.ID    = deleteID;
			swap2.type  = SAMPLE;

			swapList.push_back(swap1);
			swapList.push_back(swap2);

		}


		if ( dom.getSamplingMode()==1 || dom.getSamplingMode()==3 || dom.getSamplingMode()==5 ){ 

			if ( currentMCStep >= dom.getEquilibrationStep() ){	

				long insertID = dom.getBiasedId(SAMPLE); 
				long deleteID = dom.getRandomId(TARGET);	

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);

			}
			else{

				long insertID = dom.getRandomId(SAMPLE); 
				long deleteID = dom.getRandomId(TARGET);	

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);



			}

		}

		if ( dom.getSamplingMode()==2 || dom.getSamplingMode()==4 || dom.getSamplingMode()==6 ){ 
		
			if ( dom.getLocalBlock()==1 && dom.getBiasedBlock()==0 ){	

				long deleteID = dom.getRandomId(TARGET);	
				long insertID = dom.getLocalId(deleteID); 

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);
		
			}		
	
			if ( dom.getLocalBlock()==0 && dom.getBiasedBlock()==0 ){	

				long insertID = dom.getRandomId(SAMPLE); //choose Sample
				long deleteID = dom.getRandomId(TARGET); //choose Target

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);
			}		


			if ( dom.getLocalBlock()==0 && dom.getBiasedBlock()==1 ){ 

				long insertID = dom.getBiasedId(SAMPLE); 
				long deleteID = dom.getRandomId(TARGET);	

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);
			}
		}

		if ( dom.getSamplingMode()==7 ){ //Local moves

			if ( currentMCStep >= dom.getEquilibrationStep() ){	

				long deleteID = dom.getRandomId(TARGET);	
				long insertID = dom.getLocalId(deleteID); 

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);

			}
			else{

				long insertID = dom.getRandomId(SAMPLE); 
				long deleteID = dom.getRandomId(TARGET);	

				Swap swap1,swap2;

				swap1.ID    = insertID;
				swap1.type  = TARGET;
				swap2.ID    = deleteID;
				swap2.type  = SAMPLE;

				swapList.push_back(swap1);
				swapList.push_back(swap2);
			}

		}

		if ( dom.getSamplingMode()==8 ){ //Mono	
		
			long deleteID = dom.getRandomId(TARGET);	
			long insertID = dom.getMonoId(deleteID); 

			Swap swap1,swap2;

			swap1.ID    = insertID;
			swap1.type  = TARGET;
			swap2.ID    = deleteID;
			swap2.type  = SAMPLE;

			swapList.push_back(swap1);
			swapList.push_back(swap2);

		}



	
	}
	else{
		std::cerr<<"ERROR: SWAPLIST CONSTRUCTION NOT SUPPORTED FOR MORE THAN TWO PARTICLE"<<endl;
		exit(1);
		for(int i=0;i<getNtrialmove();i++){
			int chosenType  = getRandomNo();
			Swap swap;
			if(chosenType == TARGET)
			    swap.type = SAMPLE;
			else
				swap.type = TARGET;

			swap.ID  = dom.getRandomId(chosenType);
			swapList.push_back(swap);
		}
	}
    return swapList;
}
