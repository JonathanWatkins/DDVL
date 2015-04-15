//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryCustom.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <sstream>

//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/ini_parser.hpp>

#include "GeometryCustom.hpp"
#include "CSimulation.hpp"
#include "CParticle.hpp"
#include "FileOutput.hpp"
#include "BinnedAccumulator.hpp"


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryCustom::GeometryCustom(CSimulation * sim_)
{
	
	sim=sim_;
	
	fout = sim->GetFileOutput();
	
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
    
    omega=0;
	Amp=0;
        
}

GeometryCustom::~GeometryCustom()
{
	delete triangulatedParticlesList;
    delete triangulatedLinesList;
    delete AParticlesList;
    delete OtherParticlesList;
	
	delete Vxofy;
	delete Vxofyt;
	
    	
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
	
    //
    sim->ReadVariableFromBatchFile(Amp,"Wall.Amp");
    sim->ReadVariableFromBatchFile(omega,"Wall.omega");
    
	
	std::cout << "   Job Header loaded.\n\n";
	
	
}

void GeometryCustom::InitialiseGeometry()
{
	InitialiseFiles();
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
	
	delx = xhi - xlo;		   
	dely = yhi - ylo;
	
	Vxofy = new BinnedAccumulator(ylo,yhi,binsize);
	Vxofyt = new BinnedAccumulator(ylo,yhi,binsize);
	
	
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
		
		if (x < xlo ||	x > xhi || y < ylo || y > yhi)
		{		
			std::stringstream oss;
			oss << "GeometryCustom::CheckEscapedVortices() Vortices have escaped from the simulation box.";
			oss << " (x,y) = (" << x << ", " << y << ")";
			oss << " (x,y)_(t-1) = (" << p->get_lastx() << ", " << p->get_lasty() << ")";
			throw std::runtime_error(oss.str());
			
			//double xval = xl+(xh-xlo)*(rand() % 1000)/1000.0;
			//double yval = xl+(xh-xlo)*(rand() % 1000)/1000.0;
			//p->set_pos(xval,yval);
			
		}	
	
	}
	
	
	
	
}

void GeometryCustom::KeepParticlesInSimBoxX()
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p!=AParticlesList->end(); ++p)
	{
		
		TestX(p);
		
	}
	
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p!=OtherParticlesList->end(); ++p)
	{
		
		TestX(p);
		
	}
	
}

void GeometryCustom::KeepParticlesInSimBoxY()
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p!=AParticlesList->end(); ++p)
	{
		
		TestY(p);
		
	}
	
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p!=OtherParticlesList->end(); ++p)
	{
		
		TestY(p);		
	}
	
}

void GeometryCustom::TestX(std::list<CParticle>::iterator p)
{
	double x = p->get_x();
				
	if (x < xlo && x > -forcerange )  // left of sim box in -forcerange < x < xlo
	{	
		double testx = x+delx;
		if (testx < xlo || testx > xhi) throw std::runtime_error("GeometryCustom::KeepParticlesInSimBoxXY() Particles escaped the simulation box.");
		p->set_x(testx);
		//std::cout << testx << std::endl;
	}
			
	if (x > xhi && x < xhi+forcerange )  // right of sim box in xhi < x < xhi
	{	
		double testx = x-delx;
		if (testx < xlo || testx > xhi) throw std::runtime_error("GeometryCustom::KeepParticlesInSimBoxXY() Particles escaped the simulation box.");
		p->set_x(testx);
		//std::cout << testx << std::endl;
	}
	
	
}

void GeometryCustom::TestY(std::list<CParticle>::iterator p)
{
	double y = p->get_y();
		
	if (y < ylo && y > -forcerange )  // left of sim box in -forcerange < y < ylo
	{	
		double testy = y+dely;
		if (testy < ylo || testy > yhi) throw std::runtime_error("GeometryCustom::KeepParticlesInSimBoxXY() Particles escaped the simulation box.");
		p->set_y(testy);
	}
			
	if (y > yhi && y < yhi+forcerange )  // right of sim box in yhi < y < yhi
	{	
		double testy = y-dely;
		if (testy < ylo || testy > yhi) throw std::runtime_error("GeometryCustom::KeepParticlesInSimBoxXY() Particles escaped the simulation box.");
		p->set_y(testy);
	}

}


void GeometryCustom::KeepParticlesInSimBoxXY()
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p!=AParticlesList->end(); ++p)
	{
			TestX(p);
			TestY(p);
	}
	
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p!=OtherParticlesList->end(); ++p)
	{
			TestX(p);
			TestY(p);
	}
	
}

void GeometryCustom::AddParticlesForDT(std::list<CParticle> & iList)
{
	// Add A particles
	iList.clear();
	std::list<CParticle>::iterator it = iList.end();
	iList.insert(it,AParticlesList->begin(),AParticlesList->end());
	//return;
	iList.insert(it,OtherParticlesList->begin(),OtherParticlesList->end());
	
	if (wrapx==true && wrapy==false) WrapVorticesX(iList);
	if (wrapx==false && wrapy==true) WrapVorticesY(iList);
	if (wrapx==true && wrapy==true) WrapVorticesXY(iList);
	
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	// Temp code for testing the wrapping in the tube geomery
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	/*static bool doneoutput = false;
	
	if (doneoutput==true) return;
	
	fout->AddFileStream("wraptest","wraptest.txt");
	std::stringstream oss;
	for (std::list<CParticle>::iterator p = iList.begin();
		p != iList.end(); ++p)
	{
		oss << p->get_type() << " " << p->get_x() << " " << p->get_y() << std::endl;
		
	}
	
	doneoutput=true;
	
	fout->RegisterOutput("wraptest", oss.str());
	*/
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	// End of Temp code for testing the wrapping in the tube geomery
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
}

void GeometryCustom::PerStepAnalysis()
{
	  OutputParticlePositions(); 
	  //OutputParticleCount();
	  CalculateVxofyProfile();
	  CalculateVxofytProfile();
	  OutputVxofyEvolveProfile();
}

void GeometryCustom::EndofSimAnalysis()
{
	OutputFinalParticlePositions();
	OutputAverages();
	OutputVxofyProfile();

}

void GeometryCustom::PerStepUpdates()
{
	UserUpdates();
	if (wrapx==true && wrapy==false) KeepParticlesInSimBoxX();
	if (wrapx==false && wrapy==true) KeepParticlesInSimBoxY();
	if (wrapx==true && wrapy==true) KeepParticlesInSimBoxXY();
	CheckEscapedVortices();
}

void GeometryCustom::WrapVorticesX(std::list<CParticle>& jList)
{
		
	// Add periodic x particles	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
		DoWrapX(p, jList);
	}
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
		p!=OtherParticlesList->end(); ++p )
	{
		DoWrapX(p, jList);
	}
}

void GeometryCustom::WrapVorticesY(std::list<CParticle>& jList)
{
		
	// Add periodic y particles	
	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
		DoWrapY(p, jList);
	}
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
		p!=OtherParticlesList->end(); ++p )
	{
		DoWrapY(p, jList);
	}

}

void GeometryCustom::WrapVorticesXY(std::list<CParticle>& jList)
{
		
	// Add periodic y particles	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
			
			DoWrapXY(p, jList);			
			
	}
	
	// Add periodic y particles	
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
		p!=OtherParticlesList->end(); ++p )
	{
			
			DoWrapXY(p, jList);			
			
	}
	

}

void GeometryCustom::DoWrapX(std::list<CParticle>::iterator p, std::list<CParticle>& jList)
{
	if (p->get_x() <= xlo+forcerange)  // forcerange
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y());
		newVortex.set_type('B');
		newVortex.set_ghost();
		jList.push_back(newVortex);
	}
	if (p->get_x() >= xhi-forcerange) //channelWidth-forceRange
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y());
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}
}


void GeometryCustom::DoWrapY(std::list<CParticle>::iterator p, std::list<CParticle>& jList)
{
	if (p->get_y() <= ylo+forcerange)  // forcerange
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}
	if (p->get_y() >= yhi-forcerange) //channelWidth-forceRange
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}
}

void GeometryCustom::DoWrapXY(std::list<CParticle>::iterator p, std::list<CParticle>& jList)
{

	if (p->get_y() <= ylo+forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_y() >= yhi-forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_x() <= xlo+forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y());
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_x() >= xhi-forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y());
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	// corners
	if (p->get_y() <= ylo+forcerange && p->get_x() <=xlo+forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y()+dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_y() >= yhi-forcerange && p->get_x() >= xhi-forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y()-dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_y() <= ylo+forcerange && p->get_x() >= xhi-forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()-delx,newVortex.get_y()+dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	if (p->get_y() >= yhi-forcerange && p->get_x() <= xlo+forcerange)
	{
		CParticle newVortex;
		newVortex = (*p);
		newVortex.set_pos(newVortex.get_x()+delx,newVortex.get_y()-dely);
		newVortex.set_ghost();
		newVortex.set_type('B');
		jList.push_back(newVortex);
	}

	
}

void GeometryCustom::OutputFinalParticlePositions()
{
	
	int t = sim->get_t();
	int simulation_time = sim->get_simulation_time();
	// At the end of the simulation, output vortex positions
	if (t!=sim->get_simulation_time()+1) return;
	
	std::stringstream oss;
		
	for(std::list<CParticle>::iterator p = AParticlesList->begin();
		p != AParticlesList->end(); ++p)
	{
		oss << p->get_type() << " " << p->get_x() << " " << p->get_y();
		
		if ( std::distance(p,AParticlesList->end()) != 1 )
		{
			oss << std::endl;
		}
	}
	
	for(std::list<CParticle>::iterator p = OtherParticlesList->begin();
		p != OtherParticlesList->end(); ++p)
	{
		oss << p->get_type() << " " << p->get_x() << " " << p->get_y();
		
		if ( std::distance(p,OtherParticlesList->end()) != 1 )
		{
			oss << std::endl;
		}
	}
	
	fout->RegisterOutput("posfile", oss.str());
	
	std::cout << "Writing final vortex positions...done" << std::endl;
	
}

void GeometryCustom::OutputParticlePositions()
{
	std::stringstream oss;
	
	
	int t = sim->get_t();
	static bool header = false;	
	if (header==false)
	{
		header=true;
		fout->RegisterOutput("guifile","# This file contains frame data\n # { t, numofparticles, {id1,type1,ghost1,x1,y1,velx1,vely1,coordnum1},...,{idN, typoN, ghostN, xN,yN,velxN,velyN,coordnumN}}\n");   
	}
	
	if (t%sim->get_triangulationInterval()!=0 || t%sim->get_framedataInterval()!=0) return;
	
	// counts number of active particles

	
	int activeParticleCount=0;
	bool first=true;
	for (std::list<CParticle>::iterator p = triangulatedParticlesList->begin();
			p != triangulatedParticlesList->end(); ++p)
	{
		if (p->get_ghost()==true) continue;		
		activeParticleCount++;
		if (first==false)
		{  
			oss << ", ";
		}

		first = false;
		oss << "{"
						 << p->get_id() << ", "
						 << p->get_type() << ", "
						 << p->get_ghost() << ", "
						 << p->get_x() << ", " 
						 << p->get_y() << ", " 
						 << p->get_velx() << ", "
						 << p->get_vely() << ", "
						 << p->get_coord_num()
						 << "}";											

			
	}
	
	//std::cout << activeParticleCount << std::endl;
	
	/*for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
			p != OtherParticlesList->end(); ++p)
	{
				
		if (first==false)
		{  
			oss << ", ";
		}

		first = false;
			
		
		oss << "{"
						 << p->get_id() << ", "
						 << p->get_type() << ", "
						 << p->get_ghost() << ", "
						 << p->get_x() << ", " 
						 << p->get_y() << ", " 
						 << p->get_velx() << ", "
						 << p->get_vely()
						 << "}";											
		
			
	}*/
	
	
	oss << "}" << std::endl;;
	
	std::stringstream oss2;
	oss2 << "{" << t << ", " << activeParticleCount << ", ";
	
	oss2 << oss.str();
	
	//std::cout << oss2.str();
	fout->RegisterOutput("guifile",oss2.str()); 
}

void GeometryCustom::OutputAverages()
{	
	int t = sim->get_t();
	if (sim->get_simulation_time()+1!=t) throw std::runtime_error("Averages must be output at the end of the simulation.");
		std::stringstream oss;
		
		
		oss << "Time and space averaged quantities" << std::endl
			<< "  M2 (just stochastic term): " << sim->get_M2Average() << std::endl
			<< "  M2 (all terms): " << sim->get_M2FullAverage() << std::endl;
		fout->RegisterOutput("avfile",oss.str());
		
	std::cout << "Writing final averages...done" << std::endl;
	
}

std::list<CParticle> * GeometryCustom::GetIParticles()
{
	return AParticlesList;
}

void GeometryCustom::GetJParticles(std::list<CParticle> & jList)
{
	// Add A particles
	jList.clear();
	std::list<CParticle>::iterator it = jList.end();
	jList.insert(it,AParticlesList->begin(),AParticlesList->end());


	// Add CE particles
	jList.insert(it,OtherParticlesList->begin(),OtherParticlesList->end());

	if (wrapx==true && wrapy==false) WrapVorticesX(jList);
	if (wrapx==false && wrapy==true) WrapVorticesY(jList);
	if (wrapx==true && wrapy==true) WrapVorticesXY(jList);
	
}

void GeometryCustom::GetPeriodicity()
{
	std::string pstr;
	sim->ReadVariableFromBatchFile(pstr, "Geometry.periodicity");
	// if instr x, wrap x. If instr y, wrap y
	//std::size_t found = pstr.find("x");
	if (pstr.find("x")!=std::string::npos) wrapx=true;
	if (pstr.find("y")!=std::string::npos) wrapy=true;
} 

void GeometryCustom::OutputParticleCount()
{
	
	std::cout << "Langevin particles: " << AParticlesList->size() << std::endl;
 }
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	UserUpdates
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
 void GeometryCustom::UserUpdates()
 {
	 // Add functions here to be run every timestep
	 OscillateTopCE();
 }
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
//
//  OscillateTopCE()
//
//   V = V0 * Cos(omega * t)
//   x = x0 + V0/omega*Sin(omega*t)
//   x0 = 0
//   --> x = V0/omega*Sin(omega*t)
//   
//   xmax = V0/omega  --> if max amplitude = a0/4, V0 = 0.006 and a0 = 1
//   --> omega = V0/amp
//   phase space ->  amp = a0/10  to a0/4
//                   omega = 0.004  to 0.06
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
 
void GeometryCustom::OscillateTopCE()
{
	double t = sim->get_time();
	
	static double V0 = omega*Amp;
	
	for (std::list<CParticle>::iterator p = OtherParticlesList->begin();
		p!=OtherParticlesList->end(); ++p)
	{
		if (p->get_y() > 10)
		{	
			p->set_vel(V0*cos(omega*t), p->get_vely());
			p->set_x(p->get_x()+p->get_velx()*dt);	
		}	
	}
} 

 
 
 
void GeometryCustom::InitialiseFiles()
 {
	// add files to outputter
		
	fout->AddFileStream("posfile", "posdata.txt");
	fout->AddFileStream("guifile", "guidata.dat");
	//fout->AddFileStream("framevel", "framevel.txt");
	//fout->AddFileStream("Nd", "Nd.txt");
	fout->AddFileStream("avfile", "averagesdata.txt");
	fout->AddFileStream("Vxofy","Vxofyprofile.txt");
	fout->AddFileStream("periods","periods.txt");
	fout->AddFileStream("Vxofy_evolve","Vxofyprofile_evolve.txt");
	fout->AddFileStream("Vxofyt","Vxofytprofile.txt");
 
 }
 
void GeometryCustom::CalculateVxofyProfile()
{
	int t = sim->get_t();
	
	static int numslices = 10;
	static double period_timesteps = 2*pi/omega/dt;
	static double slice_timesteps = period_timesteps/numslices;
	
	static bool FIRSTTIME = true;
	
	std::stringstream oss;
	
	if (FIRSTTIME==true)
	{
		FIRSTTIME=false;
		oss << "period_timesteps  slice_timesteps" << std::endl;
		oss << period_timesteps << " " << slice_timesteps << std::endl;

		std::cout << "period_timesteps  slice_timesteps" << std::endl;
		std::cout << period_timesteps << " " << slice_timesteps << std::endl;
		
		fout->RegisterOutput("periods",oss.str());
		oss.str("");
	}
	
		
	//std::cout << period_timesteps << " " << slice_timesteps << std::endl;
	int sliceindex = (t % int(period_timesteps))/slice_timesteps;
	//std::cout << sliceindex << std::endl;
	oss << "t: " << t << " " << sliceindex << std::endl;
	fout->RegisterOutput("periods",oss.str());
    
		
	if (sliceindex!=1) return;
 	//std::cout << "out";
	for (std::list<CParticle>::iterator p = triangulatedParticlesList->begin();
 			p != triangulatedParticlesList->end(); ++p)
 	{
 		if (p->get_ghost()==true) continue;		
 		
 		
 		double y = p->get_y(); 
 		double f = p->get_velx();
		Vxofy->AddValue(y,f);
			
 	}
	
	
}

void GeometryCustom::CalculateVxofytProfile()
{
	int t = sim->get_t();
	
	static int numslices = 20;
	static double period_timesteps = 2*pi/omega/dt;
	static double slice_timesteps = period_timesteps/numslices;
	
	static bool FIRSTTIME = true;
	
	std::stringstream oss;
	
	if (FIRSTTIME==true)
	{
		FIRSTTIME=false;
		oss << "period_timesteps  slice_timesteps" << std::endl;
		oss << period_timesteps << " " << slice_timesteps << std::endl;
		
		std::cout << "Vxoft parameters" << std::endl;
		std::cout << "period_timesteps  slice_timesteps" << std::endl;
		std::cout << period_timesteps << " " << slice_timesteps << std::endl;
		
		fout->RegisterOutput("periods",oss.str());
		oss.str("");
	}
	
		
	//std::cout << period_timesteps << " " << slice_timesteps << std::endl;
	
	static int lastsliceindex = 0;
	
	int sliceindex = (t % int(period_timesteps))/slice_timesteps;
	if (sliceindex != lastsliceindex)
	{
		OutputVxofytProfile();
		Vxofyt->ClearValues();
		lastsliceindex = sliceindex;
	}
	
	//std::cout << sliceindex << std::endl;
	//oss << "t: " << t << " " << sliceindex << std::endl;
	//fout->RegisterOutput("periods",oss.str());
    
		
	//if (sliceindex!=1) return;
 	//std::cout << "out";
	for (std::list<CParticle>::iterator p = triangulatedParticlesList->begin();
 			p != triangulatedParticlesList->end(); ++p)
 	{
 		if (p->get_ghost()==true) continue;		
 		
 		
 		double y = p->get_y(); 
 		double f = p->get_velx();
		Vxofyt->AddValue(y,f);
			
 	}
	
	
}

void GeometryCustom::OutputVxofyProfile()
{
	std::stringstream oss;
	Vxofy->GetBinnedAverages(oss);
	fout->RegisterOutput("Vxofy",oss.str());
		
}

void GeometryCustom::OutputVxofyEvolveProfile()
{
	
	int t = sim->get_t();
	if (t%1000 == 0)	
	{
		std::cout << "VxofyEvolve output." << std::endl;
		std::stringstream oss;
		oss << t << std::endl;
		Vxofy->GetBinnedAverages(oss);
		fout->RegisterOutput("Vxofy_evolve", oss.str());	
	}
}

void GeometryCustom::OutputVxofytProfile()
{
	
	int t = sim->get_t();

	std::cout << "Vxofyt slice output." << std::endl;
	std::stringstream oss;
	oss << t << std::endl;
	Vxofyt->GetBinnedAverages(oss);
	fout->RegisterOutput("Vxofyt", oss.str());	
	
}


 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
