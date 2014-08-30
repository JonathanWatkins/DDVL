#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

// Class header
#include "CSimulation.hpp"

// Custom classes
#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CDelTriangle.hpp"
#include "CCoord.hpp"
#include "CBin.hpp"
#include "CPin.hpp"
#include "CLineIDs.hpp"
#include "CVersion.hpp"
#include "CCell.hpp"
#include "CRunningStats.hpp"
#include "delaunay.hpp"
#include "Utilities.hpp"

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
#include <boost/ptr_container/ptr_list.hpp>
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

CSimulation::CSimulation()
{
	startTime=clock();
	seedtime = time(0);
	lasttime=time(0);
	srand ( seedtime );
	running = false;
	initialised=false;
	simulation_time=1;
	simulation_initialised = false;
	paused=false;
	M2=0;
	M2Full=0;
	M2Sum=0;
	M2FullSum=0;
	avXVel=0;
	avYVel=0;
	bathLength=0;
	bathWidth=0;
	channelLength=0;
	channelWidth=0;
	sourceBfield=0;
	sinkBfield=0;
	framedataInterval=5;
	original_T=0;
	rhov=1;
	source_rhov=1;
	sink_rhov=1;	
	cellSize=0;
	binsize=0;
	Nv=0;
	Nmis=0;
	thermostat="";
	lorentzForce=0;
	jobtag="";
	triangulateReadRun=true;
	calcTrajectories=false;
	Ap=1;
	DTcount=0;
	fcount=0;
	DTtime=0;
	ftime=0;
	zoom=1;
	flat_channel_ends=false;
	reflected_channel_ends=false;
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
	lastchangedSource=0;
	lastchangedSink=0;
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
	
	calculateFinishTime();
	
	geom->ReplaceEscapedVortices();	
		
	//std::cout << "b";
	clock_t startclock = clock();
		
	integrator->Integrate();
	
	ftime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	fcount++;
	
	startclock = clock();
	delaunayTriangulation(vorticesList);
	//delVortexList=vorticesList;	
	DTtime+=(clock()-startclock)/(double)CLOCKS_PER_SEC;
	DTcount++;
			
	calculateAvVel();
	updateBathDensities();
	
}

int CSimulation::Initialise(std::string jobBatchFileLocation_)
{
	
	jobBatchFileLocation=jobBatchFileLocation_;
	
	iniread_JobBatchFile();
	
	// make jobnum from seedtime and jobtag from jobbatch file	
	std::stringstream oss;
	oss << "job" << seedtime << "-" << jobtag; 
	
	jobnum=oss.str();
	std::cout << "JobNum: " << jobnum << std::endl << std::endl;
		
	// run initial functions
	
	initialise_files();
	
	geom = CreateGeometry();
			
	geom->InitialiseVortices();
	
	geom->InitialisePins();
	
	geom->InitialiseDisorder();
	
	// set some variables
	
	vortexSize=0.2*get_a0();
	
	A=2*kB*temp/eta;
	
	// initialise integrator
	
	integrator = new CParallelEulerIntegrator(this);
	
	std::cout << "firstPin " << geom->GetFirstPin().get_x() << ", " << geom->GetFirstPin().get_y() << std::endl;
		
	iniwrite_jobheader();
	
	std::cout << "   Simulation initialised.\n\n";
		
	initialised = true;
	
	return 0;
}

void CSimulation::initialise_files()
{
	
	// make new directory for data
	std::cout << "Initialising files..." << std::endl;
	fileOutputter.setJobDirectory(jobnum);
	
	// add files to outputter
	
	//fileOutputter.addFileStream("forceterms", "forceterms.txt");
	
	if (outputType==0)  // none
	{
		std::cout << "   No results files...\n\n";
		return;
	}
	
	if (outputType==1) // positions and velocity, final positions
	{
		fileOutputter.addFileStream("guiheader", "jobheader.ini");
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

void CSimulation::iniwrite_jobheader()
{
	if (outputType==0)
		return;
	
	*(fileOutputter.getFS("guiheader")) << "[Overview]\n"
	<< "jobnum=" << jobnum << std::endl
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") <<"[ReadableBatchOptions]\n"
	<< "runtype=0" << std::endl
	<< "geometry=";
	if(geometry==channel) *fileOutputter.getFS("guiheader") << "channel";
	else if(geometry==tube) *fileOutputter.getFS("guiheader") << "tube";
	
	
	*fileOutputter.getFS("guiheader") << std::endl;
	if(geometry==channel || geometry==tube || geometry==wedge || geometry==BSCCO)
	{
		*fileOutputter.getFS("guiheader")<< "sourceBfield="<< sourceBfield << std::endl
		<< "sinkBfield="<< sinkBfield << std::endl
		<< "bathLength="<< bathLength/a0 << std::endl
		<< "bathWidth="<< bathWidth/b0 << std::endl;
	}
	
	*fileOutputter.getFS("guiheader") << "channelLength="<< channelLength/a0 << std::endl
	<< "channelWidth="<< channelWidth/b0 << std::endl
	<< "simulationTime="<< simulation_time << std::endl
	<< "temp=" << original_T << std::endl
	<< "lorentzForce=" << lorentzForce << std::endl
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") <<"[InputData]\n"
	<< "altPosFile=" << alt_pos_file << std::endl
	<< "altPosFileName="<< pos_file_name << std::endl
	<< std::endl;
	

	*fileOutputter.getFS("guiheader") <<"[BatchOptions]\n"
	<< "runtype=0" << std::endl
	<< "geometry="<< geometry << std::endl;
	if(geometry==channel || geometry==tube)
	{
		*fileOutputter.getFS("guiheader")<< "sourceBfield="<< sourceBfield << std::endl
		<< "sinkBfield="<< sinkBfield << std::endl
		<< "bathLength="<< bathLength << std::endl
		<< "bathWidth="<< bathWidth << std::endl;
	}
	
	*fileOutputter.getFS("guiheader") << "channelLength="<< channelLength << std::endl
	<< "channelWidth="<< channelWidth << std::endl
	<< "simulationTime="<< simulation_time << std::endl
	<< "temp=" << original_T << std::endl
	<< "LorentzForce=" << lorentzForce << std::endl
	
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") << "[Job]\n"
	<< "jobtag=" << jobtag << std::endl
	<< std::endl;


	*fileOutputter.getFS("guiheader") <<"[ConfigVariables]\n"
	<< "f0=" << f0 << std::endl
	
	<< std::endl;
		
	*fileOutputter.getFS("guiheader") << "[DrawingVariables]\n"
	<< "firstPinx=" << geom->GetFirstPin().get_x() << std::endl
	<< "firstPiny=" << geom->GetFirstPin().get_y() << std::endl
	<< std::endl;

	*fileOutputter.getFS("guiheader") << "[GeneralParameters]\n"
	<< "a0=" << a0 << std::endl
	<< "binSize=" << binsize << std::endl
	<< "cellSize=" << cellSize << std::endl
	<< "pi=" << pi << std::endl 
	<< "Phi=" << Phi << std::endl 
	<< "forceRange=" << forceRange << std::endl 
	<< "eta=" << eta << std::endl 
	<< "kB=" << kB << std::endl 
	<< "mu0=" << mu0 << std::endl 
	<< "lambda=" << lambda << std::endl 
	<< "Ap=" << Ap << std::endl 
	<< "dt=" << dt << std::endl 
	<< "tau=" << tau << std::endl 
	<< "drawInterval=" << framedataInterval << std::endl 
	<< "triangulationInterval=" << triangulationInterval << std::endl
	<< "thermostat=" << thermostat << std::endl
	<< "disorderDensity=" << disorderDensity << std::endl
	<< "disorderStrength=" << disorderStrength << std::endl
	<< "disorderRange=" << disorderRange << std::endl
	<< "vvForce=" << vvForce << std::endl
	<< "Av=" << Av << std::endl
	<< "Rv=" << Rv << std::endl
	<< "epsilon=" << epsilon << std::endl
	<< "sigma=" << sigma << std::endl
	
	
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") << "[BathParameters]\n"
	<< "applyBathVelocities="<< applyBathVelocities << std::endl
	<< "applyStiffBath="<< applyStiffBath << std::endl
	<< "flatChannelEnds="<< flat_channel_ends << std::endl
	<< "reflectedChannelEnds="<< flat_channel_ends << std::endl
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") << "[WallParameters]\n"
	<< "applyBounceBack="<< applyBounceBack << std::endl
	<< std::endl;
	
	*fileOutputter.getFS("guiheader") << "[DDVL]\n"
	<< "version=" << version.get_versionStr() << std::endl
	<< std::endl;
 
  
}

void CSimulation::delaunayTriangulation( std::list<CParticle> vorticesList_)
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
  
	/*if (geometry==channel) 
	{
		while (p != delLinesList.end())
		{
			bool removed=false;
		
			
			if ((p->get_x1() < bathLength && p->get_y1() < 0  &&  p->get_x2() > bathLength)
			|| (p->get_x2() < bathLength && p->get_y2() < 0  &&  p->get_x1() > bathLength))
			{
				p=delLinesList.erase(p);
				removed=true;	
			} 
			
			if ((p->get_x1() < bathLength && p->get_y1() > channelWidth+b0  &&  p->get_x2() > bathLength)
			|| (p->get_x2() < bathLength && p->get_y2() > channelWidth+b0  &&  p->get_x1() > bathLength)) 
			{
				p=delLinesList.erase(p);
				removed=true;	
			} 
			
			if ((p->get_x1() > bathLength+channelLength && p->get_y1() < 0  &&  p->get_x2() < bathLength+channelLength)
			|| (p->get_x2() > bathLength+channelLength && p->get_y2() < 0  &&  p->get_x1() < bathLength+channelLength)) 
			{
				p=delLinesList.erase(p);
				removed=true;	
			} 
			
			if ((p->get_x1() > bathLength+channelLength && p->get_y1() > channelWidth+b0  &&  p->get_x2() < bathLength+channelLength)
			|| (p->get_x2() > bathLength+channelLength && p->get_y2() > channelWidth+b0  &&  p->get_x1() < bathLength+channelLength)) 
			{
				p=delLinesList.erase(p);
				removed=true;	
			} 
			
		
			if (removed==false) { ++p; }
		}
	}
	else if (geometry==tube)
	{
		while (p != delLinesList.end())
		{
			bool removed=false;
		
			
			if (p->get_y1() <= removetopchannely || p->get_y2() <= removetopchannely ||
				p->get_y1() >=removebottomchannely || p->get_y2() >= removebottomchannely)
			{
				p=delLinesList.erase(p);
				removed=true;	
			} 
		
			if (removed==false) { ++p; }
		}
	}*/
	
}

void CSimulation::iniread_JobBatchFile() 
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
	calcTrajectories=pt.get<bool>("GeneralParameters.calcTrajectories");
	
	
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
	if(vvForce==GaussianType)
	{
		Av=pt.get<double>("Interactions.Av");
		Rv=pt.get<double>("Interactions.Rv");
	}
	else if (vvForce==BesselType)
	{
		Phi=pt.get<double>("Interactions.Phi");
		//mu0=pt.get<double>("Interactions.mu0");
		lambda=pt.get<double>("Interactions.lambda");	
	}
	else if (vvForce==LJType)
	{
		epsilon=pt.get<double>("Interactions.epsilon");
		sigma=pt.get<double>("Interactions.sigma");
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
	
	runtype=pt.get<int>("Header.runtype");
	outputType=pt.get<int>("Header.outputType");
	
	channelLength=pt.get<double>("Header.channelLength")*a0;
	channelWidth=pt.get<double>("Header.channelWidth")*b0;
		
	
	if (vvForce==BesselType || vvForce==BessLogType)
	{	
		if (geometry==channel)
		{
			sourceBfield=pt.get<double>("Header.sourceBfield");
			sinkBfield=pt.get<double>("Header.sinkBfield");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=pt.get<double>("Header.bathWidth")*b0;
		}
		else if (geometry==tube)
		{
			sourceBfield=pt.get<double>("Header.sourceBfield");
			sinkBfield=pt.get<double>("Header.sinkBfield");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=channelWidth;
		}
		else if (geometry==periodic)
		{
			Bfield=pt.get<double>("Header.Bfield");
		}
		else
		{
			std::ostringstream oss;
			oss.str("");
			oss <<"iniread_JobBatchFile() Cannot run with geometry: " << geometry << " and vvForce: " << vvForce << std::endl;
			throw std::runtime_error(oss.str());
		}
	}
	else if (vvForce==GaussianType)
	{
		if (geometry==channel)
		{
			sourceBfield=pt.get<double>("Header.sourcerhov");
			sinkBfield=pt.get<double>("Header.sinkrhov");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=pt.get<double>("Header.bathWidth")*b0;
		}
		else if (geometry==tube)
		{
			sourceBfield=pt.get<double>("Header.sourcerhov");
			sinkBfield=pt.get<double>("Header.sinkrhov");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=channelWidth;
		}
		else
		{
			std::ostringstream oss;
			oss.str("");
			oss <<"iniread_JobBatchFile() Cannot run with geometry: " << geometry << " and vvForce: " << vvForce << std::endl;
			throw std::runtime_error(oss.str());
		}
	}
	else if (vvForce==LJType)
	{
		if (geometry==channel)
		{
			source_rhov=pt.get<double>("Header.sourcerhov");
			sink_rhov=pt.get<double>("Header.sinkrhov");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=pt.get<double>("Header.bathWidth")*b0;
		}
		else if (geometry==tube)
		{
			source_rhov=pt.get<double>("Header.sourcerhov");
			sink_rhov=pt.get<double>("Header.sinkrhov");
			bathLength=pt.get<double>("Header.bathLength")*a0;
			bathWidth=channelWidth;
		}
		else
		{
			std::ostringstream oss;
			oss.str("");
			oss <<"iniread_JobBatchFile() Cannot run with geometry: " << geometry << " and vvForce: " << vvForce << std::endl;
			throw std::runtime_error(oss.str());
		}
	}
	std::cout << "   sourceBfield: " << sourceBfield << std::endl; 
	std::cout << "   sinkBfield: " << sinkBfield << std::endl; 
	
	std::cout << "   source_rhov: " << source_rhov << std::endl; 
	std::cout << "   sink_rhov: " << sink_rhov << std::endl; 
	
		
	simulation_time=pt.get<double>("Header.simulationTime");
	lorentzForce=pt.get<double>("Header.lorentzForce");  
  
  	original_T=pt.get<double>("Header.temp");  
	
	jobtag=pt.get<std::string>("Job.jobtag");  
	
	std::cout << "   Job Header loaded.\n\n";
	
}

void CSimulation::calculateFinishTime()
{
	
	if (0==t%1000)
	{
		double newtime=time(0);
		MonitorPeriod=newtime-lasttime;
		lasttime=newtime;

		double hours = floor(MonitorPeriod*(simulation_time-t)/1000.0/60.0/60.0);
		double minutes = (int)(MonitorPeriod*(simulation_time-t)/1000.0/60.0)%60;
		std::ostringstream oss;
		oss << "Estimated finish time : " << hours << "h " << minutes << "m";
		finishTimeStr = oss.str();
	
		std::cout << finishTimeStr << std::endl;
	}

	
	
}

void CSimulation::calculateAvVel()
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
		if (p->get_x()>=bathLength && p->get_x() <=bathLength+channelLength)
		{
			count++;
			spaceSumX=spaceSumX+p->get_velx();
			spaceSumY=spaceSumY+p->get_vely();
			
			
		}
		
	}
	
	
	avXVel=avXVel+spaceSumX/(double)count;
	avYVel=avYVel+spaceSumY/(double)count;
	if (outputType==1 || outputType==2) *fileOutputter.getFS("framevel") << t << " " << spaceSumX/double(count) << " " << spaceSumY/double(count) << " " << get_tAvSAvVelX() << " " << get_tAvSAvVelY() << std::endl;
	
	
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

double CSimulation::calcSinkB()
{
	
	normaliseSinkStr = "";
	
	double aaverage=0;
	int numa=0;
	for (std::list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if ( (midx > bathLength+channelLength-binsize/2.0 && midx<bathLength+channelLength+binsize/2.0) &&
				(midy>0 && midy<channelWidth)) {
			//continue;
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
											
			aaverage=aaverage+linelength;
			numa++;	
		}				
	}
	aaverage=aaverage/(double)numa;
	
	double Beff=2*Phi/(sqrt((double)3)*aaverage*aaverage);		
	std::ostringstream oss;
	oss.str("");
	
	oss << "Br(x="<<(bathLength+channelLength)/a0<<"a0)=" << Beff << "T";
	
	normaliseSinkStr = oss.str();	
	
	return Beff;
	
}

double CSimulation::calcSourceB()
{
	normaliseSourceStr = "";
	double aaverage=0;
	int numa=0;
	for (std::list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
		
		double ycutoffset = (geometry==wedge) ? channelWidth/2-(binsize/2.0*tan(pi/6)+2*b0) : channelWidth/2;
		double ymin=channelWidth/2-ycutoffset;
		double ymax=channelWidth/2+ycutoffset;				
		if ( midx> bathLength-binsize/2.0 && midx < bathLength+binsize/2.0 &&
				(midy>ymin && midy<ymax))
		{
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
											
			aaverage=aaverage+linelength;
			numa++;	
		}			
	}
	aaverage=aaverage/(double)numa;
	
	double Beff=2*Phi/(sqrt((double)3)*aaverage*aaverage);		
	std::ostringstream oss;
	oss.str("");
	
	oss << "Bl(x="<<bathLength/a0<<"a0)=" << Beff << "T";
	
	normaliseSourceStr = oss.str();	
	
	return Beff;
	
}

void CSimulation::calcBfield()
{
	BfieldStr = "";
	double aaverage=0;
	int numa=0;
	for (std::list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if ( midx > 0 && midx < channelLength &&
				(midy > 0 && midy < channelWidth))
		{
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
											
			aaverage=aaverage+linelength;
			numa++;	
		}			
	}
	aaverage=aaverage/(double)numa;
	
	double Beff=2*Phi/(sqrt((double)3)*aaverage*aaverage);		
	std::ostringstream oss;
	oss.str("");
	
	oss << "B = " << Beff << "T";
	
	BfieldStr = oss.str();	
	
	//return Beff;
	
}

double CSimulation::calcSourceRhoV()
{
	int vortex_count=0;
	
	for (std::list<CParticle>::iterator p = vorticesList.begin();
			p!=vorticesList.end(); ++p)
	{
			if (p->get_x() < bathLength)
				vortex_count++;
	}
	
	normaliseSourceStr = "";
	std::ostringstream oss;
	oss.str("");
	double unitsofarea= bathLength*bathWidth/a0/b0;
	oss << "rhov_l(x="<<bathLength/a0<<"a0)=" << double(vortex_count)/unitsofarea;
	
	normaliseSourceStr = oss.str();	
	
	return double(vortex_count)/unitsofarea;
	
}

double CSimulation::calcSinkRhoV()
{
	int vortex_count=0;
	
	for (std::list<CParticle>::iterator p = vorticesList.begin();
			p!=vorticesList.end(); ++p)
	{
			if (p->get_x() > channelLength+bathLength)
				vortex_count++;
	}
	
	normaliseSinkStr = "";
	std::ostringstream oss;
	oss.str("");
	double unitsofarea= bathLength*bathWidth/a0/b0;
	oss << "rhov_l(x="<<(bathLength+channelLength)/a0<< "a0)="<< double(vortex_count)/unitsofarea;
	
	normaliseSinkStr = oss.str();	
	
	return double(vortex_count)/unitsofarea;
	
}

CParticle CSimulation::findClosestParticle (CParticle a_)
{
  
	int j = 0;
	
	double smallest=100*a0;
	
	CParticle* result;
	
   for (std::list<CParticle>::iterator p = vorticesList.begin();
			p != vorticesList.end(); ++p)
	{
			
		double ds2 = (p->get_x() - a_.get_x())*(p->get_x() - a_.get_x())+
									(p->get_y() - a_.get_y())*(p->get_y() - a_.get_y());
      
    if (smallest > ds2)
    {
			smallest = ds2;
			result = &(*p);
    }
    
   }
   
   return (*result);
 

}

CParticle CSimulation::findClosestPin (CParticle a_)
{
	int j = 0;
	
	double smallest=1000*a0;
	
	CParticle result;
	
  for (std::list<CParticle>::iterator p = pinsList.begin();
			p != pinsList.end(); ++p)
	{
			
		double ds2 = (p->get_x() - a_.get_x())*(p->get_x() - a_.get_x())+
									(p->get_y() - a_.get_y())*(p->get_y() - a_.get_y());
      
    if (smallest > ds2)
    {
			smallest = ds2;
			result = *p;
			std::cout << smallest << std::endl;
    }
    
   }
   
   return result;
 

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

void CSimulation::updateBathDensities()
{
	if (geometry==periodic) return;  // do not run for periodic system
		
	double actualSource=0;
	double actualSink=0;
	
	double expectedSource=0;
	double expectedSink=0;
	
	
	if (vvForce==BesselType || vvForce==BessLogType)
	{
		actualSource = calcSourceB();
		actualSink = calcSinkB();
		
		expectedSource = sourceBfield;
		expectedSink = sinkBfield;
	}
	else if (vvForce==GaussianType || vvForce==LJType)
	{
		actualSource = calcSourceRhoV();
		actualSink = calcSinkRhoV();
	
		expectedSource = source_rhov;
		expectedSink = sink_rhov;
		
	}
	
	// calculate zone densities
	int sourceCount=0;
	int sinkCount=0;
	
	// relaxation_time
	double relaxation_time = a0/fabs(get_tAvSAvVelX())/channelWidth*b0/2.0;
	
	if (t%100==0) std::cout << "relax time: " << relaxation_time << std::endl;
	
	// count densities
	for (std::list<CParticle>::iterator p = vorticesList.begin();
			p != vorticesList.end(); ++p)
	{
		if (p->get_x() < bathLength && (p->get_y()<0 || p->get_y()>channelWidth))
			sourceCount++;
		
		if ( p->get_x() > bathLength+channelLength && (p->get_y()<0 || p->get_y()>channelWidth))
			sinkCount++;
	}
	
	if (dt*t-dt*lastchangedSource>=relaxation_time)		// update source after relaxation time
	{
		if (actualSource<expectedSource)
		{
			if (AddParticleToBath("source")) lastchangedSource=t;
			std::cout << "Relaxation time: " << relaxation_time << std::endl;
		}
		if (actualSource>expectedSource)
		{
			if (RemoveParticleFromBath("source")) lastchangedSource=t;
			std::cout << "Relaxation time: " << relaxation_time << std::endl;
		}
	}
	
	if (dt*t-dt*lastchangedSink>=relaxation_time)		// update sink after relaxation time
	{
		if (actualSink<expectedSink)
		{
			if (AddParticleToBath("sink")) lastchangedSink=t;
			std::cout << "Relaxation time: " << relaxation_time << std::endl;
		}
		if (actualSink>expectedSink)
		{
			if (RemoveParticleFromBath("sink")) lastchangedSink=t;
			std::cout << "Relaxation time: " << relaxation_time << std::endl;
		}
	}
	
}

bool CSimulation::AddParticleToBath(std::string location_)
{
		if (location_.compare("source") != 0 && location_.compare("sink") != 0)
			throw std::runtime_error("AddParticleToBath() location_ is not source or sink");
		
		
		
		if (location_.compare("source") == 0)
		{
			
			
			if (geometry==channel || geometry==tube)
			{
				CParticle newVortex;
				
				double xval = (2*a0)*(rand() % 1000)/1000.0;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				xval=xval+a0/2.0;
			
				newVortex.set_pos(xval,yval);
				vorticesList.push_back(newVortex);
			}
			else throw std::runtime_error ("AddParticleToBath() not valid geometry");
			
			//output newVortex added data
				
		}
		else if (location_.compare("sink") == 0)
		{
			
			
			if (geometry==channel || geometry==tube)
			{
				CParticle newVortex;
				
				double xval = (2*a0)*(rand() % 1000)/1000.0;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				xval=xval-a0/2.0;
				
				// offset to end of system
				xval = 2*bathLength + channelLength-xval;
				
				newVortex.set_pos(xval,yval);
				vorticesList.push_back(newVortex);

			}
			else throw std::runtime_error ("AddParticleToBath() not valid geometry");
			
		
				
		}		
		
		//output newVortex added data
		if (outputType!=0 && outputType!=1 ) *fileOutputter.getFS("newVortexfile") << " " << t << " " << 1 << std::endl;
		
		std::cout << "Particle added at timestep:  " << t << " in " << location_ << std::endl;
		
		return true;
}

bool CSimulation::RemoveParticleFromBath(std::string location_)
{
		if (location_.compare("source") != 0 && location_.compare("sink") != 0)
			throw std::runtime_error("RemoveParticleFromBath() location_ is not source or sink");
		
		double removalx;
		
		int vortex_to_remove = -1;
		
		// make a vector of pointers to particles in the target bath
			
		std::vector<CParticle> targetVortices;
			std::vector<CParticle> otherVortices;
		
		if (location_.compare("source") == 0)
		{
			
			
			
			removalx = geom->GetRemovalSourceX();
				
			
			
			// calculate which vortices are in removal zone
			
			for (std::list<CParticle>::iterator p = vorticesList.begin();
				p != vorticesList.end(); ++p) {
				if (p->get_x() < removalx)
				{ 
					//sinkCount++;
					targetVortices.push_back(*p);
				}
				else otherVortices.push_back(*p);
				
			}
			
			//choose a random sink vortex to be removed
	
			if (targetVortices.size()!=0)
			{
					vortex_to_remove = rand() % targetVortices.size();
			}
		}	
		else if (location_.compare("sink") == 0)
		{
			removalx = geom->GetRemovalSinkX();
			
		
			
			
			// calculate which vortices are in removal zone
			
			for (std::list<CParticle>::iterator p = vorticesList.begin();
				p != vorticesList.end(); ++p) {
				if (p->get_x() > removalx)
				{ 
					//sinkCount++;
					targetVortices.push_back(*p);
				}
				else otherVortices.push_back(*p);
				
			}
			
			//choose a random sink vortex to be removed
			
			if (targetVortices.size()!=0)
			{
					vortex_to_remove = rand() % targetVortices.size();
			}
			
		}
		
		if (vortex_to_remove!=-1)
		{ // remove a sinkVortex
			std::vector<CParticle>::iterator p = targetVortices.begin() + vortex_to_remove;
			targetVortices.erase(p);

			//update vorticesList without the removed vortex
			vorticesList.clear();
			
			std::copy( otherVortices.begin(), otherVortices.end(), std::back_inserter( vorticesList ) );
			std::copy( targetVortices.begin(), targetVortices.end(), std::back_inserter( vorticesList ) );
			
			std::cout << "Particle removed at timestep:  " << t << " from " << location_ << std::endl;
			return true;
		}
		
	return false;
		
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

