/*
 * Domain.h
 *
 */

#ifndef MC___DOMAIN_H_
#define MC___DOMAIN_H_

#include<ctime>
#include<cmath>
#include<cstdlib>
#include<random>
#include<functional>
#include <fstream>
#include <iomanip>
#include "mpi.h"
#include "Particle.h"
#include "Cell.h"
#include "UtilityMethods.h"
#include <deque>
#include <queue>
#include <algorithm>

using namespace std;

class Domain{

   public:

	  static Domain* getDomainInstance(int domainNo,IntVector3d cpuDimension,Vector3d boxSizeX,
              Vector3d boxSizeY, Vector3d boxSizeZ,IntVector3d PBC,
			  int realTypes,int totalTypes,Vector3d  sphereLayerRadii,Vector3d zoneSize,MPI_Comm communicator,
			  double cohesiveEnergy, double temperature,long nMCSteps,long writeInterval,int masterFlag,int userSeed, int SamplingMode, long StepsInBiasedBlock, long StepsInLocalBlock, double sphereRadiusTarget, int localMovesRadius, int coveringTimes,string mcSphereFile, double cylinderRadius); 

  	  // Methods

      void                     addParticle(Particle p);
      void                     setParticle(Particle p,long index);
      void                     deleteParticle(long index);
      const vector<Particle>&  getParticles(void) const;
	  void                     setCell(Cell& cell,int index);
	  void                     deleteCell(int index);

	  void                     computeCellSpeciesCount(int cellNo);
          IntVector3d               getCellSpeciesCount(int cellNo);
	  int                      getnCells(void) const;
	  int                      getDomainNo();
	  VectorRange              getDomainRange(void);
	  IntVector3d              findParticleGlobalCoordinate(Vector3d position);
	  IntVector3d              computeCellLocalCoordinate(IntVector3d globalCoordinate);
	  Cell &                   getCell(int index);
	  Cell &                   findCellwithLocalCoordinate(IntVector3d localCoordinate);
	  void                     constructCells();
	  void                     constructParticles(long nparticles, vector<long> number,
			                      vector<int> type,vector<double> mass,
							      vector<double> position,vector<double> epot, vector<int> accCounts, vector<int> rejCounts);
	  void                     createSphereNeighbors();
	  void                     createSampleNeighbors();
	  void                     setRealTypes(int realTypes);
	  int                      getRealTypes();
	  void                     setTotalTypes(int totalTypes);
	  int                      getTotalTypes();

	  double                   getTemperature(void);
	  double                   getCohesiveEnergy(void);
	  double                   getSphereRadius(void);
	  double                   getWallThickness(void);
	  Vector3d                  getBoxSizeX();
	  Vector3d                  getBoxSizeY();
	  Vector3d                  getBoxSizeZ();
          long                      getMCWriteInterval();
	  long                      getMCTotalSteps();

	  void                     setPBC(IntVector3d  PBC);
	  void                     setSphereLayer(Vector3d  sphereLayerRadii);
	  void                     createConfiguration(const vector<Swap>& swapList,vector<Particle>& sphereParticles);
	  void                     createSphere(long particleNo, vector<Particle>& sphereParticles);
	  void                     createSphereTarget(long particleNo, vector<Particle>& sphereParticles);

	  void                     shiftSphere(vector<Particle>& sphereParticles,Vector3d position);

	  long                     getnParticlesDomain();

	  void                     updateDomain(vector<Particle>& updateList , vector<Swap>& swapList);

	  void                     filterSamplingSites(void);
	  long                     getRandomId(int targetType); //TODO enum type
	  void                     writeSphereConfiguration(string filename,vector<Particle>& sphereParticles);
	  void                     writeConfiguration(string filename);
	  void                     writeCompleteConfiguration(string filename);
	  void                     writeParameterFile(string filename, int trialType,long conflicts,long finished,long swaps,long accSwaps,long rejSwaps,double totDeltaEpot,long biasedMovesPerformed,long localMovesPerformed);
	  double	           getMinimumImageSquaredDistance(Vector3d refPos,Vector3d currPos);
	  void                     updateSampleFrequency(vector<Swap>& swapList);
	  void                     addSampleSite(long mcID);
	  void                     deleteSampleSite(long mcID);
	  void                     addTargetSite(long mcID);
	  void                     deleteTargetSite(long mcID);
	  int                      getUserSeed(void);
	  long                     getTargetListSize(void) const;
      long                     getSampleListSize(void) const;
      bool                     doEnergyZoneOverlap(Vector3d v1,Vector3d v2);
      bool                     doSpheresOverlap(Vector3d v1,Vector3d v2);
      double                   getTotalEnergy(void) const;

	void increaseSuccessOrFailureFactor(vector<Particle>& Config, bool Decision);	//BIASED SAMPLING 
	double getKappa() { return Kappa_; }
	double getKappaR() { return KappaR_; }
	void setKappa(double Kappa) { Kappa_ = Kappa; }
	void setKappaR(double KappaR) { KappaR_ = KappaR; }
	double getSphereRadiusTarget() { return sphereRadiusTarget_; }
	int getLocalMovesRadius() { return localMovesRadius_; }
	long getLocalId(long targetID);
	long getMonoId(long targetID);
	long getBiasedId(int targetType);
	long biased_distribution(void);
	double getCylinderRadius() { return cylinderRadius_; }
	void setMaxSuccess( double maxSuccess) { maxSuccess_ = maxSuccess; }
	int getMaxSuccess ( void) { return maxSuccess_; }  
	void setMaxFailure( double maxFailure) { maxFailure_ = maxFailure; }
	int getMaxFailure ( void) { return maxFailure_; }  
	void setAvgSuccess( double avgSuccess) { avgSuccess_ = avgSuccess; }
	int getAvgSuccess ( void) { return avgSuccess_; }  
	void setAvgFailure( double avgFailure) { avgFailure_ = avgFailure; }
	int getAvgFailure ( void) { return avgFailure_; }  
	void setBiasStrength ( double biasStrength ){ biasStrength_ = biasStrength; }
	double getBiasStrength ( void ){ return biasStrength_; } 
	void setBellCenter ( int bellCenter ){ bellCenter_ = bellCenter; }
	int getBellCenter ( void ){ return bellCenter_; } 	
	void setInstantMaxAcc ( int instantMaxAcc ){ instantMaxAcc_ = instantMaxAcc; }
	int getInstantMaxAcc ( void ){ return instantMaxAcc_; } 
	void setInstantAvgAcc ( int instantAvgAcc ){ instantAvgAcc_ = instantAvgAcc; }
	int getInstantAvgAcc ( void ){ return instantAvgAcc_; }  
	void setInstantAvgRej ( int instantAvgRej ){ instantAvgRej_ = instantAvgRej; }
	int getInstantAvgRej ( void ){ return instantAvgRej_; }  
	void setAvgRejectionsUnderBell ( int avgRejectionsUnderBell ){ avgRejectionsUnderBell_ = avgRejectionsUnderBell; }
	int getAvgRejectionsUnderBell ( void ){ return avgRejectionsUnderBell_; } 	
	void setActiveRegionSuccess ( double activeRegionSuccess ){ activeRegionSuccess_ = activeRegionSuccess; }
	double getActiveRegionSuccess ( void ){ return activeRegionSuccess_; }
	void setActiveRegionFailure ( double activeRegionFailure ){ activeRegionFailure_ = activeRegionFailure; }
	double getActiveRegionFailure ( void ){ return activeRegionFailure_; } 
	void setInactiveRegionSuccess ( double inactiveRegionSuccess ){ inactiveRegionSuccess_ = inactiveRegionSuccess; }
	double getInactiveRegionSuccess ( void ){ return inactiveRegionSuccess_; }
	void setInactiveRegionFailure ( double inactiveRegionFailure ){ inactiveRegionFailure_ = inactiveRegionFailure; }
	double getInactiveRegionFailure ( void ){ return inactiveRegionFailure_; } 
	void setInactiveActiveRatio ( double inactiveActiveRatio ){ inactiveActiveRatio_ = inactiveActiveRatio; }
	double getInactiveActiveRatio ( void ){ return inactiveActiveRatio_; } 
	void setActiveVacancies ( double activeVacancies ){ activeVacancies_ = activeVacancies; }
	double getActiveVacancies ( void ){ return activeVacancies_; }
	double getDetailedBalanceA ( void ){ return detailedBalanceA_; } 
	void setDetailedBalanceA ( double detailedBalanceA ){ detailedBalanceA_ = detailedBalanceA; }
	double getDetailedBalanceB ( void ){ return detailedBalanceB_; } 
	void setDetailedBalanceB ( double detailedBalanceB ){ detailedBalanceB_ = detailedBalanceB; }
	double getDetailedBalanceAB ( void ){ return detailedBalanceAB_; } 
	void setDetailedBalanceAB ( double detailedBalanceAB ){ detailedBalanceAB_ = detailedBalanceAB; }
	void calculateEquilibrationStep(void);
	long getEquilibrationStep(void){ return equilibrationStep_; }
	int getBiasedBlock( void ){ return biasedBlock_; } 
	void setBiasedBlock ( int biasedBlock ){ biasedBlock_ = biasedBlock; }
	int getLocalBlock( void ){ return localBlock_; } 
	void setLocalBlock ( int localBlock ){ localBlock_ = localBlock; }
	int getSamplingMode(void){ return samplingMode_;}
	long getStepsInBiasedBlock(void){ return stepsInBiasedBlock_; }
	long getStepsInLocalBlock(void){ return stepsInLocalBlock_; }


   protected:
	  Domain(int domainNo,IntVector3d cpuDimension,Vector3d boxSizeX,
              Vector3d boxSizeY, Vector3d boxSizeZ,IntVector3d PBC,
			  int realTypes,int totalType,Vector3d  sphereLayerRadii,Vector3d zoneSize,MPI_Comm communicator,
			  double cohesiveEnergy, double temperature, long nMCSteps,long writeInterval,int masterFlag,int userSeed,int SamplingMode,long StepsInBiasedBlock,long StepsInLocalBlock, double sphereRadiusTarget, int localMovesRadius,int coveringTimes, string mcSphereFile, double cylinderRadius); 
      virtual ~Domain();

   private:
      int                   domainNo_;
      Vector3d              domainSizeX_,domainSizeY_,domainSizeZ_;
      Vector3d              totalSizeX_,totalSizeY_,totalSizeZ_;
      Vector3d              cellSize_;
      IntVector3d           domainCellArray_;
      IntVector3d           globalCellArray_;
      IntVector3d           domainPosition_;
      VectorRange           domainRange_;
      vector<Particle>      particles_;
      vector<Cell>          cells_;
      vector<IntVector3d>   cellSpecies_;
      IntVector3d           PBC_;
      double               cohesiveEnergy_;
      double               temperature_;
      Vector3d              sphereLayerRadii_;
      int                   realTypes_;
      int                   totalTypes_;
      long                  totalMCSteps_;
      long                  writeMCinterval_;
      IntVector3d           cpuDimension_;
      Vector3d              zoneSize_;

      MPI_Comm              communicator_;

      int                   userSeed_;

      // Miscellaneous members Master/Slave Parallelization
      int                   masterFlag_;
      //vector<std::reference_wrapper<Particle>> samplingSites_;

      vector<long>          samplingSites_;
      vector<long>          targetSites_;

      std::mt19937           generator_;

      int biasedBlock_;
      int localBlock_;	
      long equilibrationStep_; 
      long samplingMode_;
      long stepsInBiasedBlock_;
      long stepsInLocalBlock_;
      double Kappa_;
      double KappaR_;
      double sphereRadiusTarget_;
      int localMovesRadius_;
      int coveringTimes_; 
      int maxSuccess_;
      int maxFailure_;
      double biasStrength_;
      int avgSuccess_;
      int avgFailure_;
      int bellCenter_;
      int instantMaxAcc_;
      int instantAvgRej_;
      int instantAvgAcc_;
      int avgRejectionsUnderBell_;
      double activeRegionSuccess_;
      double activeRegionFailure_;
      double inactiveRegionSuccess_;
      double inactiveRegionFailure_;
      double inactiveActiveRatio_;
      double activeVacancies_;
      double detailedBalanceA_;
      double detailedBalanceB_;
      double detailedBalanceAB_;
      string mcSphereFile_;
      double cylinderRadius_;

};

#endif /* MC___DOMAIN_H_ */
