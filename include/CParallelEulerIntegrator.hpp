#ifndef CPARALLELEULERINTEGRATOR_HPP
#define CPARALLELEULERINTEGRATOR_HPP

#include "IntegratorBase.hpp"

#include "CParticle.hpp"
#include "CCell.hpp"
#include "CRunningStats.hpp"
#include "CSimulation.hpp"

#include <boost/function.hpp>

#define MAXLINKEDLISTSIZE 100

class CSimulation;

class CParallelEulerIntegrator : public IntegratorBase
{

public:
		
	CParallelEulerIntegrator(CSimulation * sim_);
		 
	~CParallelEulerIntegrator();
	
	void Integrate();
	
	double GetM2Average() const;
	
	double GetM2FullAverage() const;
	
	double Getdt() const {return dt;}
		
	void Initialise();
	
	double GetLambda() const { return lambda; }
	double GetForceRange() const {return forceRange; }

	
private:
	
	void ClearLinkedLists();
	
	void CreateCellLinkedLists(
		CCell ** cll_,
		std::list<CParticle> * list_);						// Creates two cell linked lists from pinsList and lastvorticesList

	void CopyCellLinkedList(
		CCell ** cllsource_,
		CCell ** clltarget_);											// Copy cell linked list

	void vvInteration(
		CParticle * p_,
		std::list<CParticle> & cell_, 
		double (&force_)[2],
		boost::function<double (double, CParallelEulerIntegrator *)> func_
		);																				// Calculates the vortex-vortex interation
	
	double forceForm
			(double dist_, bool inbath_);						// Calculates the force for a given distance
	
	void CellLinkedListToList
			(CCell** cll_,
			std::list<CParticle> * vorticesList_);	// Copy cell linked list into a std::list
	
	void temperatureInteraction
			(double (&tempForce_)[2]);													// Calculates the force due to the thermostat				

	void ApplyMaxVelocities
			(CParticle * p_,
			 double &velx_,
			 double & vely_);												// Apply max velocities

	void CheckDouble
			(double num_,
			 const std::string & varname_,
			 const std::string & source_);								// Check if double is NAN, inf or undefined 

	void CheckDuplicatePositions
			(CCell **cell_);									  	// Check for duplicate vortex positions
	
	int what_icell(std::list<CParticle>::iterator a_) const;
	
	int what_jcell(std::list<CParticle>::iterator a_) const;
	
	int what_icell(CParticle * a_) const;
	
	int what_jcell(CParticle * a_) const;
	
	double LindemanTS() const;
	
	double AndersonTS() const;
	
	double gaussianRand() const;
		
private:	
	CSimulation* sim;
	double forceRange;
	double eta;
	double dt;
	double cellSize;
	double a0;
	bool applyMaxVelocities;
	double kB;
	double tau;
	int vvForce;
	double lorentzForce;
	double f0;
	double f0bath;
	std::string thermostat;
	double b0;
	double lambda;
	double oneoverlambda;
	double temp;
		// reference variables
	double M2;
	double M2Full;
	double M2Sum;
	double M2FullSum;
	double frame_force_d;
	double frame_force_t;
	CRunningStats av_force_d;
	CRunningStats av_force_t;
	
	CCell** cll;
	CCell** lastcll;
	
	
	std::vector<double> pot_table;
	
	double step;
	double olambda;
	double rmax;
	
	double fmax;

	
	
};


#endif
