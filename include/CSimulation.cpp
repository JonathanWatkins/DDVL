#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

// Class header
#include "CSimulation.hpp"

// Custom classes
#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CDelTriangle.hpp"
#include "CBin.hpp"
#include "CLineIDs.hpp"
#include "CVersion.hpp"
#include "CCell.hpp"
#include "CRunningStats.hpp"
#include "delaunay.hpp"
#include "rv_library.hpp"

// GeometryBase Types
#include "GeometryChannel.hpp"

//#include "CParameter.hpp"

// STL classes
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <set>
#include <vector>
#include <omp.h>
#include <limits.h>
#include <iomanip>

// Boost libraries
//#include <boost/ptr_container/ptr_list.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

bool xSort (CParticle lhs_,CParticle rhs_)
{
	return (lhs_.get_x()<rhs_.get_x());
}

double CSimulation::get_M2Average() const { return M2Sum/(double)t; }

double CSimulation::get_M2FullAverage() const {	return M2FullSum/(double)t; }

double CSimulation::get_tAvSAvVelY() const { return avYVel/(double)t; } //Returns time and space average of vely of channel vortices

double CSimulation::get_tAvSAvVelX() const { return avXVel/(double)t; } // Returns time and space average of velx of channel vortices

CParticle CSimulation::get_firstPin() const { return geom->GetFirstPin(); }

double CSimulation::get_channelLength() const {	return geom->GetChannelLength(); }

double CSimulation::get_channelWidth() const { return geom->GetChannelWidth(); }

double CSimulation::get_bathLength() const { return geom->GetBathLength(); }
	
double CSimulation::get_bathWidth() const { return geom->GetBathWidth(); }

CSimulation::CSimulation()
{
	startTime=clock();
	seedtime = time(0);
	lasttime=time(0);
	srand ( seedtime );
	running = false;
	initialised=false;
	simulation_time=1;
	paused=false;
	M2=0;
	M2Full=0;
	M2Sum=0;
	M2FullSum=0;
	avXVel=0;
	avYVel=0;
	framedataInterval=5;
	cellSize=0;
	binsize=0;
	Nv=0;
	Nmis=0;
	thermostat="";
	lorentzForce=0;
	jobtag="";

	Ap=1;
	DTcount=0;
	fcount=0;
	DTtime=0;
	ftime=0;

	applyBathVelocities=false;
	applyStiffBath=false;
	applyBounceBack=false;
	applyMaxVelocities=false;
	alt_pos_file=false;
	pos_file_name="";
	alt_pins_file=false;
	pins_file_name="";
	f0_rcut_correction=0;
	f0bath_rcut_correction=0;
	
	Av=0;
	Rv=0;
	epsilon=0;
	sigma=0;
	vvForce=0;

	frame_force_t = 0;
	frame_force_d = 0;
	
	version.set_versionStr("1.0.0");	

}

void CSimulation::Run()
{
	if (initialised==false)
	{
		throw std::runtime_error("Simulation not initialised.");
	}
	
	std::cout << "Simulation running ..." << std::endl; 
	std::cout << "-------------------------------------------------" << std::endl << std::endl;
	
	for (t=1; t<=simulation_time; ++t)
	{
		
		DoStep();
		OutputVortexPositions();
		
	}
	std::cout << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
	std::cout << "Simulation finished." << std::endl << std::endl; 
	
	OutputFinalVortexPositions();
	OutputPinsList();
	OutputSimulationTimes();

}

void CSimulation::DoStep()
{
	if (0==t%100) std::cout << "t: " << t << std::endl;
	
	CalculateFinishTime();
	
	geom->ReplaceEscapedVortices();	
		
	//std::cout << "b";
	clock_t startclock = clock();
		
	integrator->Integrate();
	
	ftime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	fcount++;
	
	startclock = clock();
	DelaunayTriangulation(vorticesList);
	//delVortexList=vorticesList;	
	DTtime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	DTcount++;
			
	CalculateAvVel();
	
	geom->UpdateBathDensities();
	
}

int CSimulation::Initialise(std::string jobBatchFileLocation_)
{
	
	jobBatchFileLocation=jobBatchFileLocation_;
	
	ReadJobBatchFile();
	
	// make jobnum from seedtime and jobtag from jobbatch file	
	std::stringstream oss;
	oss << "job" << seedtime << "-" << jobtag; 
	jobnum=oss.str();
	std::cout << "JobNum: " << jobnum << std::endl << std::endl;
		
	// run initial functions
	
	InitialiseFiles();
	
	geom = CreateGeometry();
	
	geom->InitialiseVortices();
	
	geom->InitialisePins();
	
	geom->InitialiseDisorder();
	
	// calculate parameters
	
	A=2*kB*temp/eta;
	
	// initialise integrator
	
	integrator = new CParallelEulerIntegrator(this);
	
	std::cout << "firstPin " << geom->GetFirstPin().get_x() << ", " << geom->GetFirstPin().get_y() << std::endl;
		
	CopyJobBatchFile();
	
	std::cout << "Simulation initialised.\n\n";
		
	initialised = true;
	
	return 0;
}

void CSimulation::InitialiseFiles()
{
	
	// make new directory for data
	std::cout << "Initialising files..." << std::endl;
	fileOutputter.setJobDirectory(jobnum);
	
	// add files to outputter
		
	if (outputType==0)  // none
	{
		std::cout << "   No results files...\n\n";
		return;
	}
	
	if (outputType==1) // positions and velocity, final positions
	{
		fileOutputter.addFileStream("jobheader", "jobheader.ini");
		fileOutputter.addFileStream("posfile", "posdata.txt");
		fileOutputter.addFileStream("guifile", "guidata.dat");
		fileOutputter.addFileStream("trajfile", "trajectories.txt");
		fileOutputter.addFileStream("framevel", "framevel.txt");
		fileOutputter.addFileStream("CoM", "CoM.txt");
		fileOutputter.addFileStream("Nd", "Nd.txt");
		fileOutputter.addFileStream("avfile", "averagesdata.txt");
		fileOutputter.addFileStream("pinsfile", "pinsdata.txt");
		return;
	}
	
    
	std::cout << "   Files initialised.\n\n";
}

void CSimulation::OutputSimulationTimes()
{

	endTime = clock();
	
	std::cout << "Run Time: " << (endTime-startTime)/(double)CLOCKS_PER_SEC << std::endl;
	
	std::cout << "Total DTtime: " << DTtime << std::endl;
	std::cout << "Total ftime: " << ftime << std::endl;
	
	
	std::cout << "DTtime: " << DTtime/(double)DTcount << " per iteration (" << DTcount << ")" << std::endl;
	std::cout << "ftime: " << ftime/(double)fcount << " per iteration (" << fcount << ")"<< std::endl;
	
}

void CSimulation::CopyJobBatchFile()
{
	// copy jobheader to job directory and change name to jobheader.ini
	std::ifstream f1(jobBatchFileLocation, std::fstream::binary);
  
  // this is a hack
  std::ostringstream oss;
	oss << jobnum << "//jobheader.ini";
  std::ofstream f2(oss.str().c_str(), std::fstream::trunc|std::fstream::binary);
  f2 << f1.rdbuf();
	// this is a hack
 
}

void CSimulation::DelaunayTriangulation( std::list<CParticle> vorticesList_)
{
	delVortexList.clear();
	delLinesList.clear();
	
	std::list<CLineIDs> lines;
	
	geom->AddParticlesForDT(vorticesList_);
	
	/* Define input points. */

	int numberofpoints = vorticesList_.size();
  
	point2d pointlist[100000];
	
	std::vector<CParticle> vorticesVector;
	
	std::copy( vorticesList_.begin(), vorticesList_.end(), std::back_inserter( vorticesVector ) );
  
	int count=0;
	for (std::vector<CParticle>::iterator p = vorticesVector.begin();
			p!=vorticesVector.end(); ++p )
	{
		p->set_coord_num(0);
		pointlist[count].x = p->get_x();
		pointlist[count].y = p->get_y();
		count++;
	}
	
	
	int *faces = NULL;
	
	int num_faces = delaunay2d((double*)pointlist,numberofpoints,&faces);
	
	int offset = 0;
	
	//std::cout << "lines: " << lines.size() << std::endl;
	
	std::vector<int> poly;
	for(int i = 0; i < num_faces; i++ )
		{
		poly.clear();
		
			
			int num_verts = faces[offset];
			
			offset++;
			for( int j = 0; j < num_verts; j++ )
			{
				int p0 = faces[offset + j];
				int p1 = faces[offset + (j+1) % num_verts];
				
				poly.push_back(p0);
				
				
			}
			
			offset += num_verts;
		
			
			
			for (int p=0;p<poly.size()-1;p++)
			{
				
				
				
				CLineIDs newline(poly[p],poly[p+1]);
			
				lines.push_back(newline);
				
			
			}
			
					
		
		}
		
	free(faces);
	
	// check delLinesList for dupilcates
	lines.sort();
	lines.unique();
	
	for (std::list<CLineIDs>::iterator p = lines.begin();
		p!=lines.end(); ++p)
	{
	
		CDelLine newDelLine;
		newDelLine.set_points(vorticesVector[p->id1].get_x(),vorticesVector[p->id1].get_y(),vorticesVector[p->id2].get_x(),vorticesVector[p->id2].get_y());
		delLinesList.push_back(newDelLine);	
	
		vorticesVector[p->id1].coordPlusOne();
		vorticesVector[p->id2].coordPlusOne();
	}
	
	std::copy( vorticesVector.begin(), vorticesVector.end(), std::back_inserter( delVortexList ) );
  
	// remove lines between baths and channel passing over the CE
	
	std::list<CDelLine>::iterator p = delLinesList.begin();
  
}

void CSimulation::ReadJobBatchFile() 
{
	std::cout << "Loading job batch file..." << std::endl;
	std::cout << "   from " << jobBatchFileLocation << std::endl;
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(jobBatchFileLocation, pt);
	
	// geometry variables
	a0= pt.get<double>("GeneralParameters.a0");
	b0=(std::sqrt((double)3)/2.0)*a0;
	geometry=pt.get<int>("Header.geometry");
	
	// analysis variables
	binsize=pt.get<double>("GeneralParameters.binSize");
		
	// physics variables and constants
	pi=pt.get<double>("GeneralParameters.pi");
	forceRange=pt.get<double>("GeneralParameters.forceRange");
	eta=pt.get<double>("GeneralParameters.eta");
	kB=pt.get<double>("GeneralParameters.kB");
	Ap=pt.get<double>("GeneralParameters.Ap");
	
	// simulation variables
	cellSize=pt.get<double>("GeneralParameters.cellSize");
	
	dt=pt.get<double>("GeneralParameters.dt");
	tau=pt.get<double>("GeneralParameters.tau");
	triangulationInterval=pt.get<int>("GeneralParameters.triangulationInterval");
	framedataInterval=pt.get<int>("GeneralParameters.framedataInterval");
		
	
	thermostat=pt.get<std::string>("GeneralParameters.thermostat");
	
	alt_pos_file = pt.get<bool>("InputData.altPosFile");
	if (alt_pos_file == true)
	{
			pos_file_name = pt.get<std::string>("InputData.altPosFileName");
	
	}
	
	alt_pins_file = pt.get<bool>("InputData.altPinsFile");
	if (alt_pins_file == true)
	{
			pins_file_name = pt.get<std::string>("InputData.altPinsFileName");
	
	}
				
			
	// interactions
	vvForce=pt.get<double>("Interactions.vvForce");
	if (vvForce==BesselType)
	{
		Phi=pt.get<double>("Interactions.Phi");
		//mu0=pt.get<double>("Interactions.mu0");
		lambda=pt.get<double>("Interactions.lambda");	
	}
	else if (vvForce==BessLogType)
	{
		Phi=pt.get<double>("Interactions.Phi");
	}
	
	// channnel disorder
	disorderDensity=pt.get<double>("GeneralParameters.disorderDensity");
	disorderStrength=pt.get<double>("GeneralParameters.disorderStrength");
	disorderRange=pt.get<double>("GeneralParameters.disorderRange");
	
	// Job header section
	
	outputType=pt.get<int>("Header.outputType");
	
	simulation_time=pt.get<double>("Header.simulationTime");
	lorentzForce=pt.get<double>("Header.lorentzForce");  
  
  temp=pt.get<double>("Header.temp");  
	
	jobtag=pt.get<std::string>("Job.jobtag");  
	
	std::cout << "   Job Header loaded.\n\n";
	
}

void CSimulation::CalculateFinishTime()
{
	
	if (0==t%1000)
	{
		double newtime=time(0);
		MonitorPeriod=newtime-lasttime;
		lasttime=newtime;

		double hours = floor(MonitorPeriod*(simulation_time-t)/1000.0/60.0/60.0);
		double minutes = (int)(MonitorPeriod*(simulation_time-t)/1000.0/60.0)%60;
		double seconds = (int)(MonitorPeriod*(simulation_time-t)/1000.0)%60;
		std::cout << "Estimated finish time : " << hours << "h " << minutes << "m " << seconds << "s" <<  std::endl;
		
	}
	
}

void CSimulation::CalculateAvVel()
{
	/*
	 *   calculates the space and time average of the x and y velocities of the channel vortices.
	 * 	 Works for both channel system and tube system.
	 *   channel vortices are defined as not source or sink vortices.
	 *   
	 *   Current avYVel= Sum (Over t) [Sum (All channel vortices)->vely]/num channel vortices ;
	 *   To get time and space average divide by t
	 * 
	 * 	 same for x
	 */ 
	double spaceSumX=0;
	double spaceSumY=0;
	int count=0;
	for(std::list<CParticle>::iterator p = vorticesList.begin();
	    p != vorticesList.end(); p++)
	{
		if (p->get_x()>=get_bathLength() && p->get_x() <=get_bathLength()+get_channelLength())
		{
			count++;
			spaceSumX=spaceSumX+p->get_velx();
			spaceSumY=spaceSumY+p->get_vely();
			
			
		}
		
	}
	
	
	avXVel=avXVel+spaceSumX/(double)count;
	avYVel=avYVel+spaceSumY/(double)count;
	if (outputType==1) *fileOutputter.getFS("framevel") << t << " " << spaceSumX/double(count) << " " << spaceSumY/double(count) << " " << get_tAvSAvVelX() << " " << get_tAvSAvVelY() << std::endl;
	
	
}

void CSimulation::OutputFinalVortexPositions()
{
	// At the end of the simulation, output vortex positions
	if (t==simulation_time)
	{
		
		for(std::list<CParticle>::iterator p = vorticesList.begin();
			p != vorticesList.end(); ++p)
		{
			*fileOutputter.getFS("posfile") << " " << p->get_x() << " " << p->get_y();
			
			if ( std::distance(p,vorticesList.end()) != 1 )
			{
				*fileOutputter.getFS("posfile") << std::endl;
			}
		}
		
	}
}

void CSimulation::OutputPinsList()
{
	static bool PinsOutputDone = false;
	// Output the pinsList	
	if (PinsOutputDone==true)
			return;
	
	PinsOutputDone=true;
	for (std::list<CParticle>::iterator p = pinsList.begin();
				p!=pinsList.end();++p) {
		*fileOutputter.getFS("pinsfile") << " " << p->get_x() << " " << p->get_y() << std::endl;
	
	}
	
	
}

void CSimulation::OutputVortexPositions()
{
		
	if (t==1) *fileOutputter.getFS("guifile") << "# This file contains frame data\n"
																	<< "# { t, numofparticles, {x1,y1,velx1,vely1,coordnum1},...,{xN,yN,velxN,velyN,coordnumN}}" << std::endl; 
	
	if (t%triangulationInterval!=0 || t%framedataInterval!=0) return;
	
	// counts number of active particles
	int numVortices=0;
	for (std::list<CParticle>::iterator p = delVortexList.begin();
		p != delVortexList.end(); ++p)
	{
		if (p->get_ghost()!=true)
		{
			numVortices++;
		}				 
	}
	
	
	*fileOutputter.getFS("guifile") << "{" << t << ", " << numVortices << ", ";
	
	
	bool first=true;
	for (std::list<CParticle>::iterator p = delVortexList.begin();
			p != delVortexList.end(); ++p)
	{
		
		if (p->get_ghost()==true) continue;

		
		if (first==false)
		{  
			*fileOutputter.getFS("guifile") << ", ";
		}

		first = false;
		*fileOutputter.getFS("guifile") << "{"
						 << p->get_id() << ", " 
						 << p->get_x() << ", " 
						 << p->get_y() << ", " 
						 << p->get_velx() << ", "
						 << p->get_vely() << ", "
						 << p->get_coord_num()
						 << "}";											

			
	}
	
	*fileOutputter.getFS("guifile") << "}" << std::endl;
}

GeometryBase * CSimulation::CreateGeometry()
{
    //char o_type = inp_->GetOtype();
    
    switch(geometry)
    {
        case 0:  return new GeometryChannel(*this);            break;
        //case 1:  return new GeometryTube(*this);             break;
        default:   throw std::runtime_error("CSimulation:CreateOption()  Bad character");
    }
}

