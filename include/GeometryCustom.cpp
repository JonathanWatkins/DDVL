//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryCustom.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "GeometryCustom.hpp"
#include "CSimulation.hpp"
#include "CParticle.hpp"


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryCustom::GeometryCustom(CSimulation * sim_)
{
	sim=new CSimulation;
	sim=sim_;
	
	triangulatedParticlesList = new std::list<CParticle>;
    triangulatedLinesList = new std::list<CDelLine>;
    AParticlesList = new std::list<CParticle>;
    OtherParticlesList = new std::list<CParticle>;
	
	a0 = 0;
	b0 = 0;
	pos_file_name = "";
	Phi = 0;
	forcerange=0;
	dt=0;
	
	binsize = 0;
	
	xlo=0;
    ylo=0;
    xhi=0;
    yhi=0;
    
    wrapx=false;
    wrapy=false;
    
        
}

GeometryCustom::~GeometryCustom()
{
	delete triangulatedParticlesList;
    delete triangulatedLinesList;
    delete AParticlesList;
    delete OtherParticlesList;
	
}

void GeometryCustom::LoadBatchFile()
{
	std::cout << "Loading job batch file..." << std::endl;
	
	// geometry variables
	sim->ReadVariableFromBatchFile(a0, "GeneralParameters.a0");
	b0=(std::sqrt((double)3)/2.0)*a0;
		
	sim->ReadVariableFromBatchFile(Phi, "Interactions.Phi");
	sim->ReadVariableFromBatchFile(forcerange, "GeneralParameters.forceRange");
	sim->ReadVariableFromBatchFile(dt, "GeneralParameters.dt");
	
	// analysis variables
	sim->ReadVariableFromBatchFile(binsize, "GeneralParameters.binSize");
	bool posfile;
	sim->ReadVariableFromBatchFile(posfile, "InputData.altPosFile");
	if (posfile == true)
	{
			sim->ReadVariableFromBatchFile(pos_file_name, "InputData.altPosFileName");
	
	}
		
	std::cout << "   Job Header loaded.\n\n";
	
	
}

void GeometryCustom::InitialiseGeometry()
{
	LoadBatchFile();
	GetPeriodicity();
	InitialiseVortices();
	InitialiseParameters();
	
}

void GeometryCustom::InitialiseParameters()
{
	
	// calculate system parameters
		
	std::cout << "a0: " << a0
			  << "Phi: " << Phi
			  << "xlo: " << xlo
			  << "xhi: " << xhi
			  << "ylo: " << ylo
			  << "yhi: " << yhi
			  << std::endl;
			   
	 
	std::cout << "Custom geometry selected." << std::endl;
	
}

void GeometryCustom::InitialiseVortices()
{

	std::cout << "Initialising Vortices..." << std::endl;
	
	
	std::cout << "   " << pos_file_name << std::endl;
	
	std::ifstream myfile (pos_file_name.c_str());
	
	if (myfile.is_open()) // Get all particle positions from file
	{
		std::cout << "   " << "Initial Vortex Positions From File" << std::endl;
		
		
		if (myfile.good())  //get header info
		{	
			myfile >> xlo >> xhi >> ylo >> yhi;
		}
		if (xhi<=xlo) throw std::runtime_error("GeometryCustom()::InitialiseVortices() Requires xhi>xlo");
		if (yhi<=ylo) throw std::runtime_error("GeometryCustom()::InitialiseVortices() Requires yhi>ylo");
		
		double xval;
		double yval;
		char type;
		
		while ( myfile.good() )
		{
			myfile >> type >> xval >> yval;
						
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			newVortex.set_type(type);
			if (type=='A') AParticlesList->push_back(newVortex);  
			else OtherParticlesList->push_back(newVortex);
				
	
		}
		myfile.close();
	
	}
	else // Make random mobile particles and CE particles
	{	
		std::stringstream oss;
		oss << "GeometryCustom::InitialiseVortices() Could not load file " << pos_file_name;
		throw std::runtime_error(oss.str());
	}
	
	std::cout 	<< "   " << "initialiseVortices() created " << AParticlesList->size() << " Langevin vortices" 
				<< " and " << OtherParticlesList->size() << " other votices." << std::endl << std::endl;
}
                     
void GeometryCustom::CheckEscapedVortices()
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p!=AParticlesList->end(); ++p)
	{
		
		double x = p->get_x();
		double y = p->get_y();
		
		if (x <= xlo ||	x >= xhi || y <= ylo || y >= yhi)
		{		
			throw std::runtime_error("GeometryCustom::CheckEscapedVortices() Vortices have escaped from the simulation box.");
			
			//double xval = xl+(xh-xlo)*(rand() % 1000)/1000.0;
			//double yval = xl+(xh-xlo)*(rand() % 1000)/1000.0;
			//p->set_pos(xval,yval);
			
		}	
	
	}
	
	
	
	
}

void GeometryCustom::AddParticlesForDT(std::list<CParticle> & iList)
{
	// Add A particles
	iList.clear();
	std::list<CParticle>::iterator it = iList.end();
	iList.insert(it,AParticlesList->begin(),AParticlesList->end());
	
	if (wrapx==true && wrapy==false) WrapVorticesX(iList);
	if (wrapx==false && wrapy==true) WrapVorticesY(iList);
	if (wrapx==true && wrapy==true) WrapVorticesXY(iList);
	
	/*static bool once = false;
	if (once == true ) return;
	for(std::list<CParticle>::iterator p = iList.begin();
			p != iList.end(); ++p)
	{
		*sim->get_FS("wraptest") << p->get_type() << " " << p->get_x() << " " << p->get_y();
		
		if ( std::distance(p,iList.end()) != 1 )
		{
			*sim->get_FS("wraptest") << std::endl;
		}
	}
	once = true;
	*/
}

void GeometryCustom::PerStepAnalysis()
{
	  OutputParticlePositions(); 
	  //OutputParticleCount();
}

void GeometryCustom::EndofSimAnalysis()
{
	OutputFinalParticlePositions();
	OutputAverages();

}

void GeometryCustom::PerStepUpdates()
{
	CheckEscapedVortices();	
	
}

void GeometryCustom::WrapVorticesX(std::list<CParticle>& iList)
{
		
	// Add periodic x particles	
	double wrapsize = forcerange;
	double delx = xhi-xlo;
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
			if (p->get_y() <= xlo+forcerange)  // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y());
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			if (p->get_y() >= xhi-forcerange) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y());
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
	}

}

void GeometryCustom::WrapVorticesY(std::list<CParticle>& iList)
{
		
	// Add periodic y particles	
	double wrapsize = forcerange;
	double dely = yhi-ylo;
	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
			if (p->get_y() <= ylo+forcerange)  // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			if (p->get_y() >= yhi-forcerange) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
	}

}

void GeometryCustom::WrapVorticesXY(std::list<CParticle>& iList)
{
		
	// Add periodic y particles	
	double wrapsize = forcerange;
	double delx = xhi-xlo;
	double dely = yhi-ylo;
	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
						
			if (p->get_y() <= ylo+forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_y() >= yhi-forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_x() <= xlo+forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y());
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_x() >= xhi-forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y());
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			// corners
			if (p->get_y() <= ylo+forcerange && p->get_x() <=xlo+forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y()+dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_y() >= yhi-forcerange && p->get_x() >= xhi-forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y()-dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_y() <= ylo+forcerange && p->get_x() >= xhi-forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y()+dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			
			if (p->get_y() >= yhi-forcerange && p->get_x() <= xlo+forcerange)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y()-dely);
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
	}

}

void GeometryCustom::OutputFinalParticlePositions()
{
	int t = sim->get_t();
	int simulation_time = sim->get_simulation_time();
	// At the end of the simulation, output vortex positions
	if (t==sim->get_simulation_time()+1)
	{
		
		for(std::list<CParticle>::iterator p = AParticlesList->begin();
			p != AParticlesList->end(); ++p)
		{
			*sim->get_FS("posfile") << p->get_type() << " " << p->get_x() << " " << p->get_y();
			
			if ( std::distance(p,AParticlesList->end()) != 1 )
			{
				*sim->get_FS("posfile") << std::endl;
			}
		}
		
		for(std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p != OtherParticlesList->end(); ++p)
		{
			*sim->get_FS("posfile") << p->get_type() << " " << p->get_x() << " " << p->get_y();
			
			if ( std::distance(p,OtherParticlesList->end()) != 1 )
			{
				*sim->get_FS("posfile") << std::endl;
			}
		}
		
		
		
		std::cout << "Writing final vortex positions...done" << std::endl;
	}
}

void GeometryCustom::OutputParticlePositions()
{
	
	int t = sim->get_t();	
	if (t==1) *sim->get_FS("guifile") << "# This file contains frame data\n"
																	<< "# { t, numofparticles, {id1,type1,ghost1,x1,y1,velx1,vely1,coordnum1},...,{idN, typoN, ghostN, xN,yN,velxN,velyN,coordnumN}}" << std::endl; 
	
	if (t%sim->get_triangulationInterval()!=0 || t%sim->get_framedataInterval()!=0) return;
	
	// counts number of active particles
	
	*sim->get_FS("guifile") << "{" << t << ", " << AParticlesList->size() + OtherParticlesList->size() << ", ";
	
	
	bool first=true;
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p != AParticlesList->end(); ++p)
	{
				
		if (first==false)
		{  
			*sim->get_FS("guifile") << ", ";
		}

		first = false;
		*sim->get_FS("guifile") << "{"
						 << p->get_id() << ", "
						 << p->get_type() << ", "
						 << p->get_ghost() << ", "
						 << p->get_x() << ", " 
						 << p->get_y() << ", " 
						 << p->get_velx() << ", "
						 << p->get_vely()
						 << "}";											

			
	}
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p != OtherParticlesList->end(); ++p)
	{
				
		if (first==false)
		{  
			*sim->get_FS("guifile") << ", ";
		}

		first = false;
		*sim->get_FS("guifile") << "{"
						 << p->get_id() << ", "
						 << p->get_type() << ", "
						 << p->get_ghost() << ", "
						 << p->get_x() << ", " 
						 << p->get_y() << ", " 
						 << p->get_velx() << ", "
						 << p->get_vely()
						 << "}";											

			
	}
	
	
	*sim->get_FS("guifile") << "}" << std::endl;
}

void GeometryCustom::OutputAverages()
{	
	int t = sim->get_t();
	if (sim->get_simulation_time()+1!=t) throw std::runtime_error("Averages must be output at the end of the simulation.");
	
		*sim->get_FS("avfile") << "Time and space averaged quantities" << std::endl;
		*sim->get_FS("avfile") << "  M2 (just stochastic term): " << sim->get_M2Average() << std::endl;
		*sim->get_FS("avfile") << "  M2 (all terms): " << sim->get_M2FullAverage() << std::endl;
	
	std::cout << "Writing final averages...done" << std::endl;
	
}

std::list<CParticle> * GeometryCustom::GetIParticles()
{
	return AParticlesList;
}

void GeometryCustom::GetJParticles(std::list<CParticle> & iList)
{
	// Add A particles
	iList.clear();
	std::list<CParticle>::iterator it = iList.end();
	iList.insert(it,AParticlesList->begin(),AParticlesList->end());


	// Add CE particles
	iList.insert(it,OtherParticlesList->begin(),OtherParticlesList->end());

	if (wrapx==true && wrapy==false) WrapVorticesX(iList);
	if (wrapx==false && wrapy==true) WrapVorticesY(iList);
	if (wrapx==true && wrapy==true) WrapVorticesXY(iList);
	
}

void GeometryCustom::GetPeriodicity()
{
	std::string pstr;
	sim->ReadVariableFromBatchFile(pstr, "Geometry.periodicity");
	// if instr x, wrap x. If instr y, wrap y
	std::size_t found = pstr.find("x");
	if (pstr.find("x")!=std::string::npos) wrapx=true;
	if (pstr.find("y")!=std::string::npos) wrapy=true;
} 

void GeometryCustom::OutputParticleCount()
{
	
	std::cout << "Langevin particles: " << AParticlesList->size() << std::endl;
 }
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
