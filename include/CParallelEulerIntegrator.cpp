#pragma warning ( disable : 239  )  // supresses warning a bug due to icc and boost compilation
#pragma warning ( disable : 809  )  // supresses warning a bug due to icc and boost compilation


#include "CParallelEulerIntegrator.hpp"

#include <list>
#include <cilk/cilk.h>
#include <string>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

#include "rv_library.hpp" 
#include "CSimulation.hpp"
#include "GeometryBase.hpp"
#include "thermostats.hpp"


// force prototypes

double BesselsForce(const double & dist_, CParallelEulerIntegrator *integrator_);


CParallelEulerIntegrator::CParallelEulerIntegrator(CSimulation *sim_)
{
	M2=0;
	M2Full=0;
	M2Sum=0;
	M2FullSum=0;
	
	// set pointers
	sim=sim_;
		
	//initialse variables to 0
	
}

void CParallelEulerIntegrator::Initialise()
{	
	//get parameters from file
	sim->ReadVariableFromBatchFile(forceRange, "GeneralParameters.forceRange");
	sim->ReadVariableFromBatchFile(eta, "GeneralParameters.eta");
	sim->ReadVariableFromBatchFile(kB, "GeneralParameters.kB");
	sim->ReadVariableFromBatchFile(cellSize, "GeneralParameters.cellSize");
	sim->ReadVariableFromBatchFile(dt, "GeneralParameters.dt");
	sim->ReadVariableFromBatchFile(tau, "GeneralParameters.tau");
	sim->ReadVariableFromBatchFile(thermostat, "GeneralParameters.thermostat");
	sim->ReadVariableFromBatchFile(vvForce, "Interactions.vvForce");
	sim->ReadVariableFromBatchFile(lambda, "Interactions.lambda");	
	sim->ReadVariableFromBatchFile(temp, "Header.temp");  
	sim->ReadVariableFromBatchFile(lorentzForce, "Header.lorentzForce");  
	applyMaxVelocities=false;
	a0=sim->Geta0();
	oneoverlambda=1/lambda;
	// cell-linked lists on heap
	
	cll = new CCell*[MAXLINKEDLISTSIZE];
	lastcll = new CCell*[MAXLINKEDLISTSIZE];
	
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		cll[i] = new CCell[MAXLINKEDLISTSIZE];
		lastcll[i] = new CCell[MAXLINKEDLISTSIZE];
		
	}
	
	
	// initialise force
	
	step = 0.0001;
	olambda = 1/1.11;
	rmax = 6.66;
	for (long i = 1; i<=int(6.66/step)+1; ++i)
	{
		double r = step*i;
	    double f = boost::math::cyl_bessel_k(1,  r*olambda);
            
	    pot_table.push_back(f);
    }    
    
	fmax = 0;

	
	
}


CParallelEulerIntegrator::~CParallelEulerIntegrator()
{
	// Clean up cll, cllp, lastcll, lastcllp
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
    delete [] cll[i];
    delete [] lastcll[i];
    
	}
	delete [] cll;
	delete [] lastcll;
	
		
};

void CParallelEulerIntegrator::Integrate()
{
	
	M2=0;
	M2Full=0;
	
	// copy the current vorticesList to the lastvorticesList
	
	std::list<CParticle> * vorticesList = sim->get_geom()->GetIParticles();
	std::list<CParticle> lastvorticesList;
	
	
	sim->get_geom()->GetJParticles(lastvorticesList);
	
	//std::cout << "vorticesList->size(): " << vorticesList->size() << std::endl;
	//std::cout << "lastvorticesList.size(): " << lastvorticesList.size() << std::endl;
	// divide the vorticesList and lastvorticesList into cells
	
	//CreateCellLinkedLists(cll, vorticesList);
	
	CreateCellLinkedLists(lastcll, &lastvorticesList);
		
	
	// loop over all cll comparing with lastcll and cllp lists
	
	
	// make vorticesVector
	std::vector<CParticle> vorticesVector;
	std::vector<CParticle>::iterator it = vorticesVector.end();
	vorticesVector.insert(it,vorticesList->begin(),vorticesList->end());
	
	
	cilk_for (long ii = 0; ii != vorticesVector.size(); ++ii)
	{
		
		CParticle * p = &vorticesVector[ii];
		
		long i = what_icell(p);
		long j = what_jcell(p);
					
		// initialise forces to be zero
		double force[2]={0,0};
		double vortexForce[2]={0,0};
		double pinsForce[2]={0,0};
		double disorderForce[2]={0,0};
		double tempForce[2]={0,0};
		
		// calculate forces due to temperature kick 
		thermostats::Andersen(temp, kB, eta, dt, tau, tempForce);
		
		// is the vortex in the bath (so will need stiff lattice adjustment
							
		// check interation between particles
		// in this and neighbouring cells
		for(int k = i-1; k<=i+1;k++)
		{
			for(int l = j-1; l<=j+1;l++)
			{
				vvInteration(p,(*lastcll[k][l].get_cellList()),vortexForce,BesselsForce);
																							
			}
		}
		
		// set raw velocity
		
		
		
		double forcep_dx = vortexForce[0]+lorentzForce;
		double forcep_dy = vortexForce[1];
		
		double forcep_tx = tempForce[0];
		double forcep_ty = tempForce[1];
		
		//*fileOutputter.getFS("forceterms") << forcep_dx << " " << forcep_dy << " " << forcep_tx << " " << forcep_ty << std::endl;
		//std::cout << forcep_dx << "  " << std::endl;
		
		p->set_force_d_t(forcep_dx, forcep_dy, forcep_tx, forcep_ty);
							
		double velx=(forcep_dx + forcep_tx)/eta;
		double vely=(forcep_dy + forcep_ty)/eta;
		
		// Apply simulation adjustments to velocity
							
		ApplyMaxVelocities(p,velx,vely);
		
		// set velocity then check if valid
		
		p->set_vel(velx,vely);
		
		CheckDouble(p->get_velx(),"velx","calculateForces()");
		CheckDouble(p->get_vely(),"vely","calculateForces()");
		
		// set last position
		p->set_lastpos(p->get_x(),p->get_y());
		
		// set new position then check if valid
		p->set_pos(p->get_x()+p->get_velx()*dt,p->get_y()+p->get_vely()*dt);
		
		CheckDouble(p->get_x(),"x","calculateForces()");
		CheckDouble(p->get_y(),"y","calculateForces()"); 
												
					
	
	}
	
	// check for duplicate positions
	//CheckDuplicatePositions(cll);
	
	// updates vortices list
	vorticesList->clear();
	
	std::list<CParticle>::iterator itv = vorticesList->end();
	vorticesList->insert(itv,vorticesVector.begin(),vorticesVector.end());
	
	//CellLinkedListToList(cll,vorticesList);
	
	ClearLinkedLists();
	
	// average forces per particle for the system
	frame_force_d = 0;
	frame_force_t = 0;

	for (std::list<CParticle>::iterator p = vorticesList->begin();
		p != vorticesList->end(); ++p )
	{
		double forcep_dx= p->get_force_dx();
		double forcep_dy= p->get_force_dy();
		double forcep_tx= p->get_force_tx();
		double forcep_ty= p->get_force_ty();
		
		frame_force_d += sqrt(forcep_dx*forcep_dx+forcep_dy*forcep_dy);
		frame_force_t += sqrt(forcep_tx*forcep_tx+forcep_ty*forcep_ty);
		
		M2 += (forcep_tx*dt/eta)*(forcep_tx*dt/eta)+(forcep_ty*dt/eta)*(forcep_ty*dt/eta);
		M2Full += (p->get_x()-p->get_lastx())*(p->get_x()-p->get_lastx())+(p->get_y()-p->get_lasty())*(p->get_y()-p->get_lasty());
		
					
		
	}			
		
	
	
	frame_force_d/=vorticesList->size();
	frame_force_t/=vorticesList->size();
	
	av_force_d.add(frame_force_d);
	av_force_t.add(frame_force_t);
	
	M2Sum+=M2/vorticesList->size()/dt;
	M2FullSum+=M2Full/vorticesList->size()/dt;	
	
	
}

//*************************************************************************************************************
// 
//	Calculates forces on particles
//
//		Uses cell linked list, and parallelised using cilk
//
//*************************************************************************************************************

double BesselsForce(const double & dist_, CParallelEulerIntegrator * integrator_)
{
	//static double lambda = integrator_->GetLambda();
	//static double rcut = integrator_->GetForceRange();
	
	//if (dist_==0)
	//{
	//	std::cout << "zero" << std::endl;
	//	force=0.0000000001*fabs(rv::MT_rand_N());
	//}
	//else
	//{
	return 0;
	//return  boost::math::cyl_bessel_k(1,  dist_/lambda);//lambda3;// - boost::math::cyl_bessel_k(1,  rcut/thislambda)/lambda3;
		
		//double io=0;
		//double ipo=0;
		//double kpo=0;
		//rv::bessik(dist_/lambda, 1, io, force, ipo, kpo);
		
		//force = 1.0/dist_ + dist_*dist_*dist_/rcut/rcut/rcut/rcut - 2*dist_/rcut/rcut;
	//}
	
}

/*double GaussianForce(const double & dist_, const bool & inbath_, CSimulation *sim_)
{
	double Av = sim_->get_Av();
	double Rv = sim_->get_Rv();
	return Av*2*dist_*exp(-(dist_*dist_)/(double)(Rv*Rv))/Rv/Rv;
					
}*/

/*double GaussianPinForce(const double & dist_, const bool & inbath_, CSimulation *sim_)
{
	double disorderRange = sim_->get_disorderRange();
	double disorderStrength = sim_->get_disorderStrength();
	return -disorderStrength*2*dist_*exp(-(dist_*dist_)/(double)(disorderRange*disorderRange))/disorderRange/disorderRange;
					
}*/

/*double LJForce(const double & dist_, const bool & inbath_, CSimulation *sim_)
{
	double sigma = sim_->get_sigma();
	double epsilon = sim_->get_epsilon();
	return 4.0*epsilon*(12*pow(sigma,12)/pow(dist_,13) - 6*pow(sigma,6)/pow(dist_,7));
}*/

//*************************************************************************************************************
// 
// Creates two cell linked lists from pinsList and lastvorticesList
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::CreateCellLinkedLists( CCell ** cll_, std::list<CParticle> * list_ )
{

	for(std::list<CParticle>::iterator q = list_->begin();
			q != list_->end(); q++)
	{
		
		cll_[what_icell(q)][what_jcell(q)].add_particle(*q); 
		//if (what_icell(*q)==2 && what_jcell(*q)==2) std::cout << q->get_x() << ", " << q->get_y() << std::endl;	
	
	}
	
}


//*************************************************************************************************************
// 
// Copy cell linked lists
//
//*************************************************************************************************************

void CParallelEulerIntegrator::CopyCellLinkedList( CCell ** cllsource_, CCell ** clltarget_)
{

	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		for(int j = 0; j < MAXLINKEDLISTSIZE; ++j)
		{
			clltarget_[i][j]=cllsource_[i][j];
				
		}
	}
	
}


//*************************************************************************************************************
// 
// Calculates the vortex-vortex interation
//
//*************************************************************************************************************

inline float lerp(float a, float b, float f)
{
    return a + f * (b - a);
}


void CParallelEulerIntegrator::vvInteration(
		CParticle * p_,
		std::list<CParticle> & cell_,  
		double (&force_)[2],
		boost::function<double (double, CParallelEulerIntegrator *)> func_
		)
{
	
	double forceSum[2]={0,0};
	
	static double forceRangesq = forceRange*forceRange;
	
	
	double px = p_->get_x();
	double py = p_->get_y();
	/*double pxqx = 0;
	double pyqy = 0;
	double rsq = 0;
	double r = 0;
	double f = 0;
	
	double oneoverr=0;
	*/
	
    double rcut = 3.33;
    double orcutsq = 1/rcut/rcut;
	
	for(std::list<CParticle>::iterator q = cell_.begin();
									q != cell_.end(); ++q)
	{
		
		// do not check interaction between particle and itself
		//if (q->get_id()==pid)
		//		continue;
		
		double pxqx = px-q->get_x();
		double pyqy = py-q->get_y();
		 
		// r is the distance between points
		// rvector is the direction from q to p
		// r hat is the unit vector pointing from q to p
		
		double rsq = pxqx*pxqx + pyqy*pyqy;
		
		// check r is a valid number
		//if (r!=r) throw std::runtime_error ("calculateForces() r is nan"); 
		
		//if (boost::math::isinf(r)) throw std::runtime_error ("calculateForces() r is inf");

		//only include vortices closer than the forceRange cutoff								
		if (rsq > forceRangesq)
				continue;
		
		double r= sqrt( rsq );
		double oneoverr=1/r;
		
		// calculate rhat
		
		if (r==0)
			continue;
			
		double rhat[2] ={pxqx*oneoverr,pyqy*oneoverr};
		
	
		//f=func_(r,this);
		//double f = boost::math::cyl_bessel_k(1,  r*oneoverlambda);
		
		
		double f;
		if (r < step) f = pot_table[0];
		else if (r > rmax) f = 0.0;
		else
		{
			int pot_lindex=static_cast<int>(floor(r/step))-1;
		    f = lerp(pot_table[pot_lindex],pot_table[pot_lindex+1],r-(pot_lindex+1)*step);
                
		}
		
		
		
		//f = 1/r;
		//forceForm(r,inbath_);
		//double rsqorcut = r*r*orcutsq;
            
		//double f = 1.0/r+rsqorcut*rsqorcut/r-2.0*rsqorcut/r;
		
		
		forceSum[0]+=f*rhat[0];
		forceSum[1]+=f*rhat[1];
						
	}
	
	force_[0]+=forceSum[0];
	force_[1]+=forceSum[1];
	
	
	
}



//*************************************************************************************************************
// 
// Calculates the force for a given distance
//
//*************************************************************************************************************

double CParallelEulerIntegrator::forceForm(double dist_, bool inbath_)
{
	double force=0;
	
	// if in bath thislambda and thisf0 should be set to the stiff values
	double thislambda = inbath_==true ? 0.2*a0 : lambda;
	double thisf0 = inbath_==true ? f0bath : f0;
	
	if (dist_==0)
	{
		std::cout << "zero" << std::endl;
		force=0.0000000001*fabs(rv::MT_rand_N());
	}
	else
	{
		force = thisf0*boost::math::cyl_bessel_k(1,  dist_/thislambda);
	}
	
	return force;
	
}


//*************************************************************************************************************
// 
// Copy cell linked list into a std::list
//
//*************************************************************************************************************

void CParallelEulerIntegrator::CellLinkedListToList(CCell** cll_,std::list<CParticle> * vorticesList_)
{
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		for(int j = 0; j < MAXLINKEDLISTSIZE; ++j)
		{
			std::list<CParticle>::iterator it = vorticesList_->begin();
			vorticesList_->insert(it,(*cll_[i][j].get_cellList()).begin(), (*cll_[i][j].get_cellList()).end());
			//std::copy( (*cll_[i][j].get_cellList()).begin(), (*cll_[i][j].get_cellList()).end(), std::back_inserter( vorticesList_ ) );
  		
		}
	}
	
}



int CParallelEulerIntegrator::what_icell(std::list<CParticle>::iterator a_) const
{
	double i= floor((a_->get_x()-sim->get_xlo()+2*cellSize)/cellSize);
	//std::cout << i << std::endl;
	if (i!=i)
	{
		std::stringstream oss;
		oss << "what_icell(): i is undefined p->get_x()= " << a_->get_x();
		throw std::runtime_error(oss.str());
	}
	if (i<0)
	{
		std::stringstream oss;
		oss << "what_icell(): (i<0) particle outside link list grid (" << a_->get_x() << ", " << a_->get_y() << ")" << std::endl; 
		throw std::runtime_error(oss.str());
	}
	if (i>MAXLINKEDLISTSIZE-1)
	{
		std::stringstream oss;
		oss << "what_icell(): (i>MAXLINKEDLISTSIZE-1) particle outside link list grid (" << a_->get_x() << ", " << a_->get_y() << ")" << std::endl; 
		throw std::runtime_error(oss.str());
	}
	
	else return i;

}

int CParallelEulerIntegrator::what_jcell(std::list<CParticle>::iterator a_) const
{
	double j= floor((a_->get_y()-sim->get_ylo()+2*cellSize)/cellSize);
	if (j!=j)
	{
		std::stringstream oss;
		oss << "what_icell(): j is undefined p->get_x()= " << a_->get_x();
		throw std::runtime_error(oss.str());
	}
	if (j<0) throw std::runtime_error("what_jcell(): (j<0) particle outside link list grid");
	
	if (j>MAXLINKEDLISTSIZE-1) throw std::runtime_error("what_jcell(): (j>MAXLINKEDLISTSIZE-1) particle outside link list grid");
	
	else return j;

}

int CParallelEulerIntegrator::what_icell(CParticle * a_) const
{
	double i= floor((a_->get_x()-sim->get_xlo()+2*cellSize)/cellSize);
	//std::cout << i << std::endl;
	if (i!=i)
	{
		std::stringstream oss;
		oss << "what_icell(): i is undefined p->get_x()= " << a_->get_x();
		throw std::runtime_error(oss.str());
	}
	if (i<0)
	{
		std::stringstream oss;
		oss << "what_icell(): (i<0) particle outside link list grid (" << a_->get_x() << ", " << a_->get_y() << ")" << std::endl; 
		throw std::runtime_error(oss.str());
	}
	if (i>MAXLINKEDLISTSIZE-1)
	{
		std::stringstream oss;
		oss << "what_icell(): (i>MAXLINKEDLISTSIZE-1) particle outside link list grid (" << a_->get_x() << ", " << a_->get_y() << ")" << std::endl; 
		throw std::runtime_error(oss.str());
	}
	
	else return i;

}

int CParallelEulerIntegrator::what_jcell(CParticle *a_) const
{
	double j= floor((a_->get_y()-sim->get_ylo()+2*cellSize)/cellSize);
	if (j!=j)
	{
		std::stringstream oss;
		oss << "what_icell(): j is undefined p->get_x()= " << a_->get_x();
		throw std::runtime_error(oss.str());
	}
	if (j<0) throw std::runtime_error("what_jcell(): (j<0) particle outside link list grid");
	
	if (j>MAXLINKEDLISTSIZE-1) throw std::runtime_error("what_jcell(): (j>MAXLINKEDLISTSIZE-1) particle outside link list grid");
	
	else return j;

}

//*************************************************************************************************************
// 
// Apply Max Velocities
//
//		If velocity is over 5b0/dt in the sink or 4b0/dt everywhere else set to maxvel
//
//*************************************************************************************************************

void CParallelEulerIntegrator::ApplyMaxVelocities(CParticle * p_, double &velx_, double & vely_)
{		
			if (applyMaxVelocities==false)
				return;
			static int num_corrections=0;
			static int num_checks=0;
			static int thist=sim->get_t();
			
			num_checks++;
			
			if (sim->get_t()==thist+2)
			{
				num_corrections=0;
				num_checks=0;
			
				//std::cout << thist << ", " << t << std::endl;
				//std::cout << "Number of vel corections: " << num_corrections << " (" << num_corrections << "/" << num_checks << ")"<< std::endl;
				thist=sim->get_t();
			}
			
			double maxvel =  0.5*a0/dt;
			 
			if(velx_>maxvel)
			{
				velx_=maxvel;
				num_corrections++;
			}
		
			if(velx_<-maxvel)
			{
				velx_=-maxvel;
				num_corrections++;
			}
		
			if(vely_>maxvel)
			{
				vely_=maxvel;
				num_corrections++;
			}
		
			if(vely_<-maxvel)
			{
				vely_=-maxvel;
				num_corrections++;
			}
			
}

//*************************************************************************************************************
// 
// Check if double is NAN, inf or undefined 
//
//		varname_ is the name of the variable being check
//		source_ is the name of the function that called CheckDouble
//
//*************************************************************************************************************


void CParallelEulerIntegrator::CheckDouble(double num_, const std::string & varname_, const std::string & source_)
{
	if (num_!=num_)
	{
		std::stringstream oss;
		oss << "CheckValidity() called from " << source_ << " : " << varname_ << " " << num_; 
		throw std::runtime_error(oss.str());
	}
	
	if (boost::math::isinf(num_))
	{
		std::stringstream oss;
		oss << "CheckValidity() called from " << source_ << " : " << varname_ << " " << num_; 
		throw std::runtime_error(oss.str());
	}		
}

//*************************************************************************************************************
// 
//	Check for duplicate vortex positions
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::CheckDuplicatePositions(CCell **cell_)
{
	for(int i = 1; i < MAXLINKEDLISTSIZE-1; ++i)
	{
		for(int j = 1; j < MAXLINKEDLISTSIZE-1; ++j)
		{
			for(std::list<CParticle>::iterator p = (*cell_[i][j].get_cellList()).begin();
								p != (*cell_[i][j].get_cellList()).end(); ++p)
			{
					for(int k = i-1; k<=i+1;k++)
					{
						for(int l = j-1; l<=j+1;l++)
						{
							for(std::list<CParticle>::iterator q = (*cell_[k][l].get_cellList()).begin();
									q != (*cell_[k][l].get_cellList()).end(); ++q)
							{
								// do not check with itself
								if (p->get_id() == q->get_id()) continue;
								
								if (p->get_x() == q->get_x() && p->get_y() == q->get_y())
								{
									p->set_x(p->get_x()+0.001*a0);
									q->set_x(q->get_x()-0.001*a0);
									std::cout << "Duplicate particle positions found and fixed" << std::endl;
								}
						
							
							}
						}
					}
				
			}
		}
	}
}

void CParallelEulerIntegrator::ClearLinkedLists()
{
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		for(int j = 0; j < MAXLINKEDLISTSIZE; ++j)
		{
			cll[i][j].clearlist();
			lastcll[i][j].clearlist();
		}
	}
	
	
}

double CParallelEulerIntegrator::GetM2Average() const {return M2Sum/sim->get_t();}
	
double CParallelEulerIntegrator::GetM2FullAverage() const {return M2FullSum/sim->get_t();}
	
