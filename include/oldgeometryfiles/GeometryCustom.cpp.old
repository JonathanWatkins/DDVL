//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryChannel.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include "GeometryChannel.hpp"

#include "CSimulation.hpp"
#include "CParticle.hpp"

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>



//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryChannel::GeometryChannel(CSimulation & sim_)
:		sim(sim_)
,   vorticesList(sim_.get_vorticesList())
,   delLinesList(sim_.get_delLinesList())
,		a0(sim_.get_a0())
,		b0(sim_.get_b0())
,		dt(sim_.get_dt())
,		binsize(sim.get_binsize())
,   	pos_file_name(sim_.GetPosFileName())
,		jobBatchFileLocation(sim_.get_jobBatchFileLocation())
{

	LoadBatchFile();
	
	
	// calculate system parameters
	xlo=sim->get_xlo();
	xhi=sim->get_xhi();
	ylo=sim->get_ylo();
	yhi=sim->get_yhi();
	
	periodicity=sim->get_periodicity(); //string containing nothing or x or y or x,y 
	// check if periodic
	
	
	
	
	
	std::cout << "Channel geometry selected." << std::endl; 
	
}

void GeometryChannel::LoadBatchFile()
{
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(jobBatchFileLocation, pt);
	
	xlo=pt.get<double>("Geometry.xlo")*a0;
	xhi=pt.get<double>("Geometry.xhi")*b0;
	ylo=pt.get<double>("Geometry.ylo")*a0;
	yhi=pt.get<double>("Geometry.yhi")*a0;
	
}

void GeometryChannel::InitialiseVortices() const
{

	std::cout << "Initialising Vortices..." << std::endl;
	std::cout << "   " << "sourceDensity: " << sourceDensity << std::endl;
	std::cout << "   " << "sinkDensity: " << sinkDensity << std::endl;
  
	std::cout << "   " << pos_file_name << std::endl;
			
	std::ifstream myfile (pos_file_name.c_str());
	
	
	if (myfile.is_open())
	{
		std::cout << "   " << "Initial Vortex Positions From File" << std::endl;
		//file = true;
		double xval;
		double yval;
		char type;
		while ( myfile.good() )
		{
			myfile >> type >> xval >> yval;
						
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			newVortex.set_type(type);
			vorticesList->push_back(newVortex);	
		}

		myfile.close();

	}
	else
	{
		std::stringstream oss;
		oss << "GeometryChannel::InitialiseVorices() Cannot load vortex positions from file " << pos_file_name << std::endl;
		throw std::runtime_error(oss.str());
	}
		
   
	std::cout << "   " << "initialiseVortices() created " << vorticesList->size() << " vortices." << std::endl << std::endl;
	
}
         
void GeometryChannel::ReplaceEscapedVortices() const
{
 	for (std::list<CParticle>::iterator p = vorticesList->begin();
			p!=vorticesList->end(); ++p)
	{
		
		if (p->get_x() <= xlo || p->get_y() <= ylo || 
				p->get_x() >= xhi  || p->get_y() >= yhi )
			std::runtime_error("Particle escaped from simulation")
	}
	
	
}

void GeometryChannel::InitialisePins()
{
	
	
	
}

void GeometryChannel::AddParticlesForDT(std::list<CParticle> & vorticesList_) const
{

	
}

void GeometryChannel::WrapSystem() const
{}

void GeometryChannel::InitialiseDisorder() const
{
	
}

double GeometryChannel::GetRemovalSourceX() const
{
 return 0;
} 

double GeometryChannel::GetRemovalSinkX() const
{
	return 0;
} 

CParticle GeometryChannel::GetFirstPin() const
{
	return firstPin;
}



void GeometryChannel::UpdateBathDensities() const
{
	
	
}

bool GeometryChannel::AddParticleToBath(std::string location_) const
{
	return false;
}

bool GeometryChannel::RemoveParticleFromBath(std::string location_) const
{
	return false;
}


double GeometryChannel::calcSinkB() const
{
	return 0;
}

double GeometryChannel::calcSourceB() const
{
	
	return 0;
}

void GeometryChannel::WrapVortices(std::list<CParticle>& vorticesList_) const
{
	return;
	
}

void GeometryCh

 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
