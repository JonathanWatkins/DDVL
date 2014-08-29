#include "CParallelEulerIntegrator.hpp"
#include "Utilities.hpp" 
#include "CSimulation.hpp"

#include <list>
#include <cilk/cilk.h>
#include <string>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

// force prototypes
double BesselsForce(double dist_, bool inbath_, CSimulation *sim_);
double BessLogForce(double dist_, bool inbath_, CSimulation *sim_);
double GaussianForce(double dist_, bool inbath_, CSimulation *sim_);
double GaussianPinForce(double dist_, bool inbath_, CSimulation *sim_);
double LJForce(double dist_, bool inbath_, CSimulation *sim_);


template<class T>
double gen_normal_3(T &generator)
{
  double rn=0;
	rn =generator();
  return rn; 
}


// Using seed for mt as time(0) does not allow runs to be repeated with same seeds
boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),
              boost::normal_distribution<>());



CParallelEulerIntegrator::CParallelEulerIntegrator(CSimulation *sim_) :
vorticesList(sim_->vorticesList), t(sim_->t), temp(sim_->temp),
M2(sim_->M2), M2Full(sim_->M2Full), M2Sum(sim_->M2Sum), 
M2FullSum(sim_->M2FullSum), frame_force_d(sim_->frame_force_d),
frame_force_t(sim_->frame_force_t), av_force_d(sim_->av_force_d), 
av_force_t(sim_->av_force_t)

{
	
	
	// set pointers
	sim=sim_;
	
	//get parameters from sim class
	applyStiffBath=sim->get_applyStiffBath();
	applyMaxVelocities=sim->get_applyMaxVelocities();
	applyBathVelocities=sim->get_applyBathVelocities();
	applyBounceBack=sim->get_applyBounceBack();
	
	geometry=sim->get_geometry();
	vvForce=sim->get_vvForce();
	
	pinsList=*sim->get_pinsList();
	disorderList=*sim->get_disorderList();
	
	dt=sim->get_dt();
	eta=sim->get_eta();
	kB=sim->get_kB();
	tau=sim->get_tau();
	Ap=sim->get_Ap();
	lorentzForce=sim->get_lorentzForce();
	
	cellSize=sim->get_cellSize();
	channelLength=sim->get_channelLength();
	channelWidth=sim->get_channelWidth();
	bathLength=sim->get_bathLength();
	bathWidth=sim->get_bathWidth();
	a0=sim->get_a0();
	flat_channel_ends=sim->get_flat_channel_ends();
	reflected_channel_ends=sim->get_reflected_channel_ends();
	f0=sim->get_f0();
	f0bath=sim->get_f0bath();
	thermostat=sim->get_thermostat();
	b0=sim->get_b0();
	
	reboundy0=sim->get_reboundy0();
	reboundy1=sim->get_reboundy1(); 
	reboundx0=sim->get_reboundx0();
	reboundx1=sim->get_reboundx1(); 
	
	bouncebacky0=sim->get_bouncebacky0();
	bouncebacky1=sim->get_bouncebacky1(); 
	bouncebackx0=sim->get_bouncebackx0();
	bouncebackx1=sim->get_bouncebackx1();
	
	lambda=sim->get_lambda();
		
	// set variables
	forceRange=sim->get_forceRange();
 	temp_f0_rcut_correction=BesselsForce(forceRange,false,sim);
	temp_f0bath_rcut_correction=applyStiffBath==true ? BesselsForce(forceRange,true,sim) : BesselsForce(forceRange,false,sim); 
	f0_rcut_correction = temp_f0_rcut_correction;
	f0bath_rcut_correction = temp_f0bath_rcut_correction;
	
	//std::cout << "   f0 rcut correction: " << f0_rcut_correction << std::endl;
	//std::cout << "   f0bath rcut correction: " << f0bath_rcut_correction << std::endl;
	
	//std::cout << "   LJ test: " << LJForce(1,false,this) << std::endl;
  
}

void CParallelEulerIntegrator::Integrate()
{
	// get this step values
	M2=0;
	M2Full=0;
	
	// copy the current vorticesList to the lastvorticesList
	std::list<CParticle> lastvorticesList=vorticesList;
	
	if (geometry==tube) Utilities::wrapVortices(lastvorticesList, sim->get_channelWidth(), sim->get_forceRange());
	if (geometry==periodic) Utilities::wrapVorticesPeriodic(lastvorticesList, sim->get_channelLength(), sim->get_channelWidth(), sim->get_forceRange());
		
	// cell-linked lists on heap
	
	CCell** cll = new CCell*[MAXLINKEDLISTSIZE];
	CCell** cllp = new CCell*[MAXLINKEDLISTSIZE];
	CCell** clldis = new CCell*[MAXLINKEDLISTSIZE];
	CCell** lastcll = new CCell*[MAXLINKEDLISTSIZE];
	CCell** nextcll = new CCell*[MAXLINKEDLISTSIZE];
	
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		cll[i] = new CCell[MAXLINKEDLISTSIZE];
		cllp[i] = new CCell[MAXLINKEDLISTSIZE];
		clldis[i] = new CCell[MAXLINKEDLISTSIZE];
		lastcll[i] = new CCell[MAXLINKEDLISTSIZE];
		nextcll[i] = new CCell[MAXLINKEDLISTSIZE];
		
	}
		
		
	// divide the vorticesList and lastvorticesList into cells
	
	CreateCellLinkedLists(cll, vorticesList);
	
	CreateCellLinkedLists(lastcll, lastvorticesList);
	
	CreateCellLinkedLists(cllp, pinsList);
	CreateCellLinkedLists(clldis, disorderList);
	
	// loop over all cll comparing with lastcll and cllp lists
	
	
	cilk_for(int i = 1; i < MAXLINKEDLISTSIZE-1; ++i)
	{
		cilk_for(int j = 1; j < MAXLINKEDLISTSIZE-1; ++j)
		{
				//if ((*cll[i][j].get_cellList()).size()!=0)
				//		std::cout << (*cll[i][j].get_cellList()).size() <<std::endl;
			
				for(std::list<CParticle>::iterator p = (*cll[i][j].get_cellList()).begin();
									p != (*cll[i][j].get_cellList()).end(); ++p)
				{
					
				
					// initialise forces to be zero
					double force[2]={0,0};
					double vortexForce[2]={0,0};
					double pinsForce[2]={0,0};
					double disorderForce[2]={0,0};
					double tempForce[2]={0,0};
					double reflectedvortexForce[2]={0,0};
					double channelEndsForce=ChannelEndsInteration(p);
					
					// calculate forces due to temperature kick 
					temperatureInteraction(p,tempForce);
					
					// initialise stress terms
					double JyyK=0;  // kinetic term
					double JyyV=0;  // potential term
		
					double JxxK=0;  // kinetic term
					double JxxV=0;  // potential term
					
					double JxyK=0;  // kinetic term
					double JxyV=0;  // potential term
					
					double JyxK=0;  // kinetic term
					double JyxV=0;  // potential term
					
					JyyK=0;// -1*q->get_vely()*q->get_vely();
					JxxK=0;// -1*q->get_velx()*q->get_velx();
					JxyK=0;// -1*q->get_velx()*q->get_vely();
					JyxK=0;// -1*q->get_velx()*q->get_vely();
								
					
					// is the vortex in the bath (so will need stiff lattice adjustment
					bool inbath=false;
					if (applyStiffBath==true)
					{
						if (p->get_x()< bathLength || p->get_x() >bathLength+channelLength)
							inbath=true;
					}				
					// check interation between particles
					// in this and neighbouring cells
					for(int k = i-1; k<=i+1;k++)
					{
						for(int l = j-1; l<=j+1;l++)
						{
							// calculate forces and stresses between vortices lastcll
							
							if (vvForce==BesselType)
							{
								vvInteration(p,(*lastcll[k][l].get_cellList()),vortexForce,JxyV,JyxV,JxxV,JyyV,inbath,BesselsForce);
								vvInteration(p,(*cllp[k][l].get_cellList()),pinsForce,JxyV,JyxV,JxxV,JyyV,inbath,BesselsForce);
								rvvInteration(p,(*lastcll[k][l].get_cellList()),reflectedvortexForce,JxyV,JyxV,JxxV,JyyV,inbath,0,BesselsForce);
							}
							else if (vvForce==BessLogType)
							{
								vvInteration(p,(*lastcll[k][l].get_cellList()),vortexForce,JxyV,JyxV,JxxV,JyyV,inbath,BessLogForce);
								vvInteration(p,(*cllp[k][l].get_cellList()),pinsForce,JxyV,JyxV,JxxV,JyyV,inbath,BessLogForce);
								rvvInteration(p,(*lastcll[k][l].get_cellList()),reflectedvortexForce,JxyV,JyxV,JxxV,JyyV,inbath,0,BessLogForce);
							}
							else if (vvForce==GaussianType)
							{
								vvInteration(p,(*lastcll[k][l].get_cellList()),vortexForce,JxyV,JyxV,JxxV,JyyV,inbath,GaussianForce);
								vvInteration(p,(*cllp[k][l].get_cellList()),pinsForce,JxyV,JyxV,JxxV,JyyV,inbath,GaussianForce);
								rvvInteration(p,(*lastcll[k][l].get_cellList()),reflectedvortexForce,JxyV,JyxV,JxxV,JyyV,inbath,0,GaussianForce);
							}
							else if (vvForce==LJType)
							{
								vvInteration(p,(*lastcll[k][l].get_cellList()),vortexForce,JxyV,JyxV,JxxV,JyyV,inbath,LJForce);
								vvInteration(p,(*cllp[k][l].get_cellList()),pinsForce,JxyV,JyxV,JxxV,JyyV,inbath,LJForce);
								rvvInteration(p,(*lastcll[k][l].get_cellList()),reflectedvortexForce,JxyV,JyxV,JxxV,JyyV,inbath,0,LJForce);
							}
							
							
							// calculate forces and stresses between quenced disorder 
							
							vvInteration(p,(*clldis[k][l].get_cellList()),disorderForce,JxyV,JyxV,JxxV,JyyV,inbath,GaussianPinForce);
												
						}
					}
					
					// set raw velocity
					
					
					
					double forcep_dx = vortexForce[0]+Ap*pinsForce[0]+disorderForce[0]+lorentzForce+channelEndsForce+reflectedvortexForce[0];
					double forcep_dy = vortexForce[1]+Ap*pinsForce[1]+disorderForce[1]+reflectedvortexForce[1];
					
					double forcep_tx = tempForce[0];
					double forcep_ty = tempForce[1];
					
					//*fileOutputter.getFS("forceterms") << forcep_dx << " " << forcep_dy << " " << forcep_tx << " " << forcep_ty << std::endl;
					
					
					p->set_force_d_t(forcep_dx, forcep_dy, forcep_tx, forcep_ty);
										
					double velx=(forcep_dx + forcep_tx)/eta;
					double vely=(forcep_dy + forcep_ty)/eta;
					
					// Apply simulation adjustments to velocity
										
					ApplyMaxVelocities(p,velx,vely);
					
					ApplyBathVelocities(p,velx,vely);			
					
					ApplyReboundConditions(p,velx,vely);
					
					ApplyBounceBackConditions(p,velx,vely);
					
					//if (velx>a0/50/dt || velx<-a0/50/dt) std::cout << "Displacement too large |dx| = " << velx*dt << " > a0/50" << std::endl;					
					//if (vely>a0/50/dt || vely<-a0/50/dt) std::cout << "Displacement too large |dy| = " << vely*dt << " > a0/50" << std::endl;					
					
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
					
					// set force 
					//p->set_force(force[0],force[1]);
					
					// calculate dx2
					//M2+=tempForce[0]*tempForce[0];
					//M2Full+=p->get_velx()*dt*p->get_velx()*dt;
					
					// Set particles local stress
					
					p->set_Jyy(JyyK+JyyV);
					p->set_Jxx(JxxK+JxxV);
					p->set_Jxy(JxyK+JxyV);
					p->set_Jyx(JyxK+JyxV);
					
					
				}
			
			
		}
	}
	
	// check for duplicate positions
	//CopyCellLinkedList(cll, nextcll);
	
	CheckDuplicatePositions(cll);
	
	// updates vortices list
	vorticesList.clear();
	CellLinkedListToList(cll,vorticesList);
	
	// Update moments
	/*if ("Anderson"==thermostat)
    {
			M2=M2/dt/Nc;
			M2Full=M2Full/dt/vorticesList.size();
		}
	else if ("Lindeman"==thermostat) M2=M2/Nc;
	else M2=0;
	*/
	//M2Sum=M2Sum+M2;
	//M2FullSum=M2FullSum+M2Full;

		
	// average forces per particle for the system
	frame_force_d = 0;
	frame_force_t = 0;

	for (std::list<CParticle>::iterator p = vorticesList.begin();
		p != vorticesList.end(); ++p )
	{
		double forcep_dx= p->get_force_dx();
		double forcep_dy= p->get_force_dy();
		double forcep_tx= p->get_force_tx();
		double forcep_ty= p->get_force_ty();
		
		frame_force_d += sqrt(forcep_dx*forcep_dx+forcep_dy*forcep_dy);
		frame_force_t += sqrt(forcep_tx*forcep_tx+forcep_ty*forcep_ty);
		
		M2 += (forcep_tx*dt/eta)*(forcep_tx*dt/eta);
		M2Full += (p->get_x()-p->get_lastx())*(p->get_x()-p->get_lastx());
		
					
		
	}			
		
	
	
	frame_force_d/=vorticesList.size();
	frame_force_t/=vorticesList.size();
	
	av_force_d.add(frame_force_d);
	av_force_t.add(frame_force_t);
	
	M2Sum+=M2/vorticesList.size()/dt;
	M2FullSum+=M2Full/vorticesList.size()/dt;	
		
		
	// Clean up cll, cllp, lastcll, lastcllp
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
    delete [] cll[i];
    delete [] cllp[i];
    delete [] lastcll[i];
    delete [] clldis[i];
    delete [] nextcll[i];
    
	}
	delete [] cll;
	delete [] cllp;
	delete [] lastcll;
	delete [] clldis;
	delete [] nextcll;
	
	
}

//*************************************************************************************************************
// 
//	Besseling paper log form of potential
//
//  	Cannot have inbath interaction, as this requires lambda=1.11
//
//*************************************************************************************************************

double BessLogForce(double dist_, bool inbath_, CSimulation *sim_)
{
	double force=0;
	
	double rcut = sim_->get_forceRange();
	
	//double thislambda = sim_->get_lambda();
	
	//double lambda2=thislambda*thislambda;
	
	if (dist_==0)
	{
		std::cout << "zero" << std::endl;
		force=0.0000000001*fabs(Utilities::gaussianRand());
	}
	//else if (dist_< 0.5*sim_->get_a0())
	//{
		//force = thisf0*boost::math::cyl_bessel_k(1,  0.5*sim_->get_a0()/thislambda)-thisf0_rcut_correction;
	//}
	else
	{
		//force = thisf0*boost::math::cyl_bessel_k(1,  dist_/thislambda)-thisf0_rcut_correction;
		double factor= (1- (dist_/rcut) * (dist_/rcut));
		
		force = (1.0/dist_) * factor * factor;
	}
	
	return force;
}


//*************************************************************************************************************
// 
//	Calculates forces on particles
//
//		Uses cell linked list, and parallelised using cilk
//
//*************************************************************************************************************

double BesselsForce(double dist_, bool inbath_, CSimulation *sim_)
{
	double force=0;
	// if in bath thislambda and thisf0 should be set to the stiff values
	double thislambda = inbath_==true ? 0.2*sim_->get_a0() : sim_->get_lambda();
	//double thisf0 = inbath_==true ? sim_->get_f0bath() : sim_->get_f0();
	//double thisf0_rcut_correction= inbath_==true ? sim_->get_f0bath_rcut_correction() : sim_->get_f0_rcut_correction();
	double rcut = sim_->get_forceRange();
	
	//double lambda3=thislambda*thislambda*thislambda;
	
	if (dist_==0)
	{
		std::cout << "zero" << std::endl;
		force=0.0000000001*fabs(Utilities::gaussianRand());
	}
	//else if (dist_< 0.5*sim_->get_a0())
	//{
	//	force = thisf0*boost::math::cyl_bessel_k(1,  0.5*sim_->get_a0()/thislambda)-thisf0_rcut_correction;
	//}
	else
	{
		force = boost::math::cyl_bessel_k(1,  dist_/thislambda);//lambda3;// - boost::math::cyl_bessel_k(1,  rcut/thislambda)/lambda3;
		//force = 1.0/dist_ + dist_*dist_*dist_/rcut/rcut/rcut/rcut - 2*dist_/rcut/rcut;
	}
	
	return force;
}

double GaussianForce(double dist_, bool inbath_, CSimulation *sim_)
{
	double Av = sim_->get_Av();
	double Rv = sim_->get_Rv();
	return Av*2*dist_*exp(-(dist_*dist_)/(double)(Rv*Rv))/Rv/Rv;
					
}

double GaussianPinForce(double dist_, bool inbath_, CSimulation *sim_)
{
	double disorderRange = sim_->get_disorderRange();
	double disorderStrength = sim_->get_disorderStrength();
	return -disorderStrength*2*dist_*exp(-(dist_*dist_)/(double)(disorderRange*disorderRange))/disorderRange/disorderRange;
					
}

double LJForce(double dist_, bool inbath_, CSimulation *sim_)
{
	double sigma = sim_->get_sigma();
	double epsilon = sim_->get_epsilon();
	return 4.0*epsilon*(12*pow(sigma,12)/pow(dist_,13) - 6*pow(sigma,6)/pow(dist_,7));
}

//*************************************************************************************************************
// 
// Creates two cell linked lists from pinsList and lastvorticesList
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::CreateCellLinkedLists( CCell ** cll_, std::list<CParticle> & list_ )
{

	for(std::list<CParticle>::iterator q = list_.begin();
			q != list_.end(); q++)
	{
		
		cll_[what_icell(*q)][what_jcell(*q)].add_particle(*q); 
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

void CParallelEulerIntegrator::vvInteration(
		std::list<CParticle>::iterator p_,
		std::list<CParticle> & cell_,  
		double (&force_)[2],
		double & JxyV_,
		double & JyxV_,
		double & JxxV_,
		double & JyyV_,
		const bool &inbath_,
		boost::function<double (double,bool, CSimulation *)> func_
		)
{
	double forceSum[2]={0,0};
	
	for(std::list<CParticle>::iterator q = cell_.begin();
									q != cell_.end(); ++q)
	{
		
		// do not check interaction between particle and itself
		if (q->get_id()==p_->get_id())
				continue;
		
		// r is the distance between points
		// rvector is the direction from q to p
		// r hat is the unit vector pointing from q to p
		
		double r= sqrt((double)
					(p_->get_x()-q->get_x())*(p_->get_x()-q->get_x())
				+ (p_->get_y()-q->get_y())*(p_->get_y()-q->get_y()));
		
		// check r is a valid number
		if (r!=r) throw std::runtime_error ("calculateForces() r is nan"); 
		
		if (boost::math::isinf(r)) throw std::runtime_error ("calculateForces() r is inf");

		//only include vortices closer than the forceRange cutoff								
		if (r > forceRange)
				continue;
		
		// calculate rhat
		double rvector[2]={p_->get_x()-q->get_x(),p_->get_y()-q->get_y()};
		double modr=sqrt((double)rvector[0]*rvector[0]+rvector[1]*rvector[1]);									
		if (modr==0)
			continue;
			
		double rhat[2] ={rvector[0]/modr,rvector[1]/modr};
		
		
		double f=func_(r,inbath_,sim);
		//forceForm(r,inbath_);
		
		forceSum[0]=forceSum[0]+f*rhat[0];
		forceSum[1]=forceSum[1]+f*rhat[1];
						
		JyyV_+=0.5*(p_->get_y()-q->get_y())*f*rhat[1];
		JxxV_+=0.5*(p_->get_x()-q->get_x())*f*rhat[0];
		JxyV_+=0.5*(p_->get_y()-q->get_y())*f*rhat[0];
		JxyV_+=0.5*(p_->get_x()-q->get_x())*f*rhat[1];
		
	}
	
	force_[0]=force_[0]+forceSum[0];
	force_[1]=force_[1]+forceSum[1];
		
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
		force=0.0000000001*fabs(Utilities::gaussianRand());
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

void CParallelEulerIntegrator::CellLinkedListToList(CCell** cll_,std::list<CParticle> & vorticesList_)
{
	for(int i = 0; i < MAXLINKEDLISTSIZE; ++i)
	{
		for(int j = 0; j < MAXLINKEDLISTSIZE; ++j)
		{
			std::copy( (*cll_[i][j].get_cellList()).begin(), (*cll_[i][j].get_cellList()).end(), std::back_inserter( vorticesList_ ) );
  		
		}
	}
	
}


int CParallelEulerIntegrator::what_icell(const CParticle & a_) const
{
	double i= floor((a_.get_x()-sim->get_firstPin().get_x()+2*cellSize)/cellSize);
	//std::cout << i << std::endl;
	std::stringstream oss;
	oss << "what_icell(): i is undefined p->get_x()= " << a_.get_x();
	if (i!=i) throw std::runtime_error(oss.str());
	oss.str("");
	if (i<0)
	{
		 oss << "what_icell(): (i<0) particle outside link list grid (" << a_.get_x() << ", " << a_.get_y() << ")" << std::endl; 
	
		 throw std::runtime_error(oss.str());
	}
	if (i>MAXLINKEDLISTSIZE-1)
	{
		 oss << "what_icell(): (i>MAXLINKEDLISTSIZE-1) particle outside link list grid (" << a_.get_x() << ", " << a_.get_y() << ")" << std::endl; 
		 throw std::runtime_error(oss.str());
	}
	else return i;

}

int CParallelEulerIntegrator::what_jcell(const CParticle & a_) const
{
	double j= floor((a_.get_y()-sim->get_firstPin().get_y()+2*cellSize)/cellSize);
	std::stringstream oss;
	oss << "what_icell(): j is undefined p->get_x()= " << a_.get_x();
	if (j!=j) throw std::runtime_error(oss.str());
	if (j<0) throw std::runtime_error("what_jcell(): (j<0) particle outside link list grid");
	if (j>MAXLINKEDLISTSIZE-1) throw std::runtime_error("what_jcell(): (j>MAXLINKEDLISTSIZE-1) particle outside link list grid");
	else return j;

}

double CParallelEulerIntegrator::AndersonTS() const
{
	double p=dt/tau;
	if (p>rand()/(double)RAND_MAX)
	{	
		
		return sqrt(2*temp*kB*eta/dt/p)*gen_normal_3(generator);
	}
	else
	{
		return 0;
			
	}
}

double CParallelEulerIntegrator::LindemanTS() const
{
	return sqrt((double)temp)*gen_normal_3(generator);
	
}


//*************************************************************************************************************
// 
// Calculates the force due to the thermostat
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::temperatureInteraction(std::list<CParticle>::iterator & q_,double (&tempForce_)[2])
{
	
	//if (q_->get_x()>bathLength && q_->get_x()<bathLength+channelLength)
		//{
		if ("Lindeman"==thermostat){
			
			tempForce_[0]=LindemanTS();
			tempForce_[1]=LindemanTS();
		}
		else if ("Anderson"==thermostat)
		{
			tempForce_[0]=AndersonTS();
			tempForce_[1]=AndersonTS();
		}
		
		if (tempForce_[0]!=tempForce_[0] || tempForce_[1] != tempForce_[1]) {
			std::cout << "t: " << t << "temp nan" << "(" << tempForce_[0] << ", " << tempForce_[1] << ")" << std::endl;
			tempForce_[0]=0;
			tempForce_[1]=0;
		}
		if (boost::math::isinf(tempForce_[0]) || boost::math::isinf(tempForce_[1]))
		std::cout << "t: " << t << "temperature inf" << "(" << tempForce_[0] << ", " << tempForce_[1] << ")" << std::endl;
		//}
}


//*************************************************************************************************************
// 
// Apply Bath Velocities
//
//		Increase sink velocity, decrease source velocity.
//		This process helps with diffussion in the sink and stops
//		vortices going too fast when added to the source
//
//*************************************************************************************************************

void CParallelEulerIntegrator::ApplyBathVelocities(std::list<CParticle>::iterator p_, double & velx_, double & vely_)
{
	if (applyBathVelocities==false)
		return; 
	if (p_->get_x() >bathLength+channelLength)
	{ // rescaled viscosity for sink vortices
		velx_=velx_*2;
		vely_=vely_*2;
	}
	else if (p_->get_x() <bathLength)
	{ // rescaled viscosity for sink vortices
		velx_=velx_*0.5;
		vely_=vely_*0.5;
	}
}


//*************************************************************************************************************
// 
// Apply Max Velocities
//
//		If velocity is over 5b0/dt in the sink or 4b0/dt everywhere else set to maxvel
//
//*************************************************************************************************************

void CParallelEulerIntegrator::ApplyMaxVelocities(std::list<CParticle>::iterator p_, double &velx_, double & vely_)
{		
			if (applyMaxVelocities==false)
				return;
			static int num_corrections=0;
			static int num_checks=0;
			static int thist=t;
			
			num_checks++;
			
			if (t==thist+2)
			{
				num_corrections=0;
				num_checks=0;
			
				//std::cout << thist << ", " << t << std::endl;
				std::cout << "Ratio of corections: " << double(num_corrections)/vorticesList.size() << " (" << num_corrections << "/" << num_checks << ")"<< std::endl;
				thist=t;
			}
			
			double maxvel = (p_->get_x()<bathLength || p_->get_x()>channelLength+bathLength) ? 0.5*b0/dt : maxvel=0.5*b0/dt;
			 
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
// Apply rebound conditions at CE
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::ApplyReboundConditions(std::list<CParticle>::iterator p_, double & velx_, double & vely_)
{
	return;
	if (geometry!=channel && geometry!=BSCCO && geometry!=wedge)
		return;
		
	if (p_->get_y()+vely_*dt>reboundy1)
	{
		double reflectedy = reboundy1-(p_->get_y()+vely_*dt-reboundy1);
		vely_=(reflectedy-p_->get_y())/dt;
	}
	else if (p_->get_y()+vely_*dt<reboundy0)
	{
		double reflectedy = reboundy0-(p_->get_y()+vely_*dt-reboundy0);
		vely_=(reflectedy-p_->get_y())/dt;
	}

	if (p_->get_x()+velx_*dt>reboundx1)
	{
		double reflectedx = reboundx1-(p_->get_x()+velx_*dt-reboundx1);
		velx_=(reflectedx-p_->get_x())/dt;
	}
	else if (p_->get_x()+velx_*dt<reboundx0)
	{
		double reflectedx = reboundx0-(p_->get_x()+velx_*dt-reboundx0);
		velx_=(reflectedx-p_->get_x())/dt;
	}



}


//*************************************************************************************************************
// 
// Apply rebound conditions at CE
//
//
//*************************************************************************************************************

void CParallelEulerIntegrator::ApplyBounceBackConditions(std::list<CParticle>::iterator p_, double & velx_, double & vely_)
{
	if (applyBounceBack==false)
		return;
		
	if (geometry!=channel && geometry!=BSCCO && geometry!=wedge)
		throw std::runtime_error("ApplyBounceBackConditions() Can only apply rebound conditions in channel, wedge or BSCCO geometry");
	
	// check if vortex hits CE. If so bounceback.
				
	if (p_->get_y()+vely_*dt>bouncebacky1)
	{
		//double reflectedy = bouncebacky1-(p_->get_y()+vely_*dt-bouncebacky1);
		//vely_=(reflectedy-p_->get_y())/dt;
		velx_=0;
		vely_=0;
	}
	else if (p_->get_y()+vely_*dt<bouncebacky0)
	{
		//double reflectedy = bouncebacky0-(p_->get_y()+vely_*dt-bouncebacky0);
		//vely_=(reflectedy-p_->get_y())/dt;
		velx_=0;
		vely_=0;
	}

	if (p_->get_x()+velx_*dt>bouncebackx1)
	{
		//double reflectedx = bouncebackx1-(p_->get_x()+velx_*dt-bouncebackx1);
		//velx_=(reflectedx-p_->get_x())/dt;
		velx_=0;
		vely_=0;
	}
	else if (p_->get_x()+velx_*dt<bouncebackx0)
	{
		//double reflectedx = bouncebackx0-(p_->get_x()+velx_*dt-bouncebackx0);
		//velx_=(reflectedx-p_->get_x())/dt;
		velx_=0;
		vely_=0;
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

void CParallelEulerIntegrator::CheckDouble(double num_, std::string varname_, std::string source_)
{
	std::stringstream oss;
	oss << "CheckValidity() called from " << source_ << " : " << varname_ << " " << num_; 
	if (num_!=num_) throw std::runtime_error(oss.str());
	if (boost::math::isinf(num_)) throw std::runtime_error(oss.str());
					
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

//*************************************************************************************************************
// 
//	ChannelEndsInteraction
//
//		Just applys flat wall potential at sink end of the channel
//
//*************************************************************************************************************

double CParallelEulerIntegrator::ChannelEndsInteration(
		std::list<CParticle>::iterator p_)
{
	if (flat_channel_ends==false) return 0;
		
	double r=p_->get_x();
	double rprime=2*bathLength+channelLength-r;
	//return 1e-9*exp(-(r*r)/1.0)/a0/a0-1e-9*exp(-(rprime*rprime)/1.0)/a0/a0;
	//return BesselsForce(r,false,this)-BesselsForce(rprime,false,this);
	return -BesselsForce(rprime,false,sim);
}

//*************************************************************************************************************
// 
// Calculates the reflected vortex-vortex interation at the source end of the channel
//
//*************************************************************************************************************

void CParallelEulerIntegrator::rvvInteration(
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
		)
{
	if (reflected_channel_ends==false) return;
	
	if (p_->get_x() > forceRange)
		return;
		
	double forceSum[2]={0,0};
	
	for(std::list<CParticle>::iterator q = cell_.begin();
									q != cell_.end(); ++q)
	{
		
		// do not check interaction between particle and itself
		if (q->get_id()==p_->get_id())
				continue;
		
		// r is the distance between points
		// rvector is the direction from q to p
		// r hat is the unit vector pointing from q to p
		
		double r= sqrt((double)
					(p_->get_x()-(wall_position_-q->get_x()))*(p_->get_x()-(wall_position_-q->get_x()))
				+ (p_->get_y()-q->get_y())*(p_->get_y()-q->get_y()));
		
		// check r is a valid number
		if (r!=r) throw std::runtime_error ("calculateForces() r is nan"); 
		
		if (boost::math::isinf(r)) throw std::runtime_error ("calculateForces() r is inf");

		//only include vortices closer than the forceRange cutoff								
		if (r > forceRange)
				continue;
		
		// calculate rhat
		double rvector[2]={p_->get_x()-(wall_position_-q->get_x()),p_->get_y()-q->get_y()};
		double modr=sqrt((double)rvector[0]*rvector[0]+rvector[1]*rvector[1]);									
		if (modr==0)
			continue;
			
		double rhat[2] ={rvector[0]/modr,rvector[1]/modr};
		
		
		double f=func_(r,inbath_,sim);
		//forceForm(r,inbath_);
		
		forceSum[0]=forceSum[0]+f*rhat[0];
		forceSum[1]=forceSum[1]+f*rhat[1];
						
		JyyV_+=0.5*(p_->get_y()-q->get_y())*f*rhat[1];
		JxxV_+=0.5*(p_->get_x()-(wall_position_-q->get_x()))*f*rhat[0];
		JxyV_+=0.5*(p_->get_y()-q->get_y())*f*rhat[0];
		JxyV_+=0.5*(p_->get_x()-(wall_position_-q->get_x()))*f*rhat[1];
		
	}
	
	force_[0]=force_[0]+forceSum[0];
	force_[1]=force_[1]+forceSum[1];
	
	
	
}


