#ifndef CPARALLELEULERINTEGRATOR_HPP
#define CPARALLELEULERINTEGRATOR_HPP

//#include "CSimulation.hpp"
#include "CParticle.hpp"
#include "CCell.hpp"
#include "CRunningStats.hpp"

#include <boost/function.hpp>


class CSimulation;

class CParallelEulerIntegrator
{

public:
		
	CParallelEulerIntegrator(CSimulation * sim_);
		 
	~CParallelEulerIntegrator();
	
	void Integrate();
	
private:
	
	void Integrate2();
	
	void ClearLinkedLists();
	
	void CreateCellLinkedLists(
		CCell ** cll_,
		std::list<CParticle> & list_);						// Creates two cell linked lists from pinsList and lastvorticesList

	void CopyCellLinkedList(
		CCell ** cllsource_,
		CCell ** clltarget_);											// Copy cell linked list

	void vvInteration(
		std::list<CParticle>::iterator p_,
		std::list<CParticle> & cell_, 
		double (&force_)[2],
		double & JxyV_,
		double & JyxV_,
		double & JxxV_,
		double & JyyV_,
		const bool &inbath_,
		boost::function<double (double,bool, CSimulation *)> func_
		);																				// Calculates the vortex-vortex interation
	
	void vvInteration2(
		CParticle *  p_,
		CParticle *  q_,
		double (&force_)[2],
		const bool &inbath_,
		boost::function<double (double,bool, CSimulation *)> func_
	);
	
	double forceForm
			(double dist_, bool inbath_);						// Calculates the force for a given distance
	
	void CellLinkedListToList
			(CCell** cll_,
			std::list<CParticle> & vorticesList_);	// Copy cell linked list into a std::list
	
	void temperatureInteraction
			(double (&tempForce_)[2]);													// Calculates the force due to the thermostat				

	void ApplyBathVelocities
			(std::list<CParticle>::iterator p_, 
			 double & velx_,
			 double & vely_);												// Apply bath velocties
	
	void ApplyMaxVelocities
			(std::list<CParticle>::iterator p_,
			 double &velx_,
			 double & vely_);												// Apply max velocities

	void ApplyMaxVelocities2(CParticle * p_, double &velx_, double & vely_);


	void CheckDouble
			(double num_,
			 std::string varname_,
			 std::string source_);								// Check if double is NAN, inf or undefined 

	void CheckDuplicatePositions
			(CCell **cell_);									  	// Check for duplicate vortex positions
	
	double ChannelEndsInteration(							
		std::list<CParticle>::iterator p_);			// If flat wall return force from channel ends
		
	void rvvInteration(		
		std::list<CParticle>::iterator p_,
		std::list<CParticle> & cell_,  
		double (&force_)[2],
		double & JxyV_,
		double & JyxV_,
		double & JxxV_,
		double & JyyV_,
		const bool &inbath_,
		double wall_position_,
		boost::function<double (double,bool, CSimulation *)> func_
		);																				// Calculates the reflected vortex-vortex interation

	int what_icell(std::list<CParticle>::iterator a_) const;
	
	int what_jcell(std::list<CParticle>::iterator a_) const;
	
	double LindemanTS() const;
	
	double AndersonTS() const;
	
	double gaussianRand() const;
	
private:	
	CSimulation* sim;
	
	double forceRange;
	
	double temp_f0_rcut_correction;
	
	double temp_f0bath_rcut_correction;
	
	bool applyStiffBath;
	
	double f0_rcut_correction;
	
	double f0bath_rcut_correction;
	
	std::list<CParticle> pinsList;
	
	std::list<CParticle> disorderList;
	
	double eta;
	
	double dt;
	
	double cellSize;
	
	double channelLength;
	
	double channelWidth;
	
	double bathLength;
	
	double bathWidth;
	
	double a0;
	
	bool applyBathVelocities;
	
	bool applyMaxVelocities;
	
	double kB;
	
	double tau;
	
	int vvForce;
	
	double Ap;
	
	double lorentzForce;
	
	double f0;
	
	double f0bath;
	
	std::string thermostat;
	
	double b0;
	
	double lambda;
	
	// reference variables
		
	std::list<CParticle> &vorticesList;
	
	double & M2;
	
	double & M2Full;
	
	double & M2Sum;
	
	double & M2FullSum;
	
	double & frame_force_d;
	
	double & frame_force_t;
	
	CRunningStats & av_force_d;
	
	CRunningStats & av_force_t;
	
	int & t;
	
	double & temp;
	
	
	CCell** cll;
	CCell** cllp;
	CCell** clldis;
	CCell** lastcll;

	
	
};


#endif
