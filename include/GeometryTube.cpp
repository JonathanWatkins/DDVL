//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryTube.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "GeometryTube.hpp"
#include "CSimulation.hpp"
#include "CParticle.hpp"


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryTube::GeometryTube(CSimulation & sim_)
:	sim(sim_)
{
	bathLength = 0;
	bathWidth = 0;
	channelLength = 0;
	channelWidth = 0;
	sourceBfield = 0;
	sinkBfield = 0;
	Phi = 0;
	a0 = 0;
	b0 = 0;
	dt = 0;
	forcerange = 0;
	pos_file_name = "";
	pins_file_name = "";
	jobBatchFileLocation = "";
	
	sourceDensity = 0;
	sinkDensity = 0;
	channelDensity = 0;
	
	removesourcex = 0;
	removesinkx = 0;
	
	removetopchannely = 0;
	removebottomchannely = 0;
	
	etchsourcex = 0;
	etchsinkx = 0;
		  
	binsize = 0;
	
	avXVel=0;
    avYVel=0;
    
    xlo=0;
    ylo=0;
    xhi=0;
    yhi=0;
        
	
}

void GeometryTube::LoadBatchFile()
{
	std::cout << "Loading job batch file..." << std::endl;
	std::cout << "   from " << jobBatchFileLocation << std::endl;
	
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(jobBatchFileLocation, pt);
	
	// geometry variables
	a0= pt.get<double>("GeneralParameters.a0");
	b0=(std::sqrt((double)3)/2.0)*a0;
		
	channelLength=pt.get<double>("Geometry.channelLength")*a0;
	channelWidth=pt.get<double>("Geometry.channelWidth")*b0;
	sourceBfield=pt.get<double>("Geometry.sourceBfield");
	sinkBfield=pt.get<double>("Geometry.sinkBfield");
	bathLength=pt.get<double>("Geometry.bathLength")*a0;
	bathWidth=pt.get<double>("Geometry.bathWidth")*b0;
	
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
	Phi=pt.get<double>("Interactions.Phi");
	//mu0=pt.get<double>("Interactions.mu0");
	lambda=pt.get<double>("Interactions.lambda");	
	
	
	// channnel disorder
	disorderDensity=pt.get<double>("GeneralParameters.disorderDensity");
	disorderStrength=pt.get<double>("GeneralParameters.disorderStrength");
	disorderRange=pt.get<double>("GeneralParameters.disorderRange");
	
	// Job header section
	
	outputType=pt.get<int>("Header.outputType");
	
	lorentzForce=pt.get<double>("Header.lorentzForce");  
  
	temp=pt.get<double>("Header.temp");  
	
	std::cout << "   Job Header loaded.\n\n";
	
	
}

void GeometryTube::InitialiseGeometry()
{
	LoadBatchFile();
	InitialiseParameters();
	InitialiseVortices();
	InitialiseCE();
}

void GeometryTube::InitialiseParameters()
{
	
	// calculate system parameters
		
	etchsourcex = 0; 
	etchsinkx = bathLength+channelLength+bathLength; 
	
	removesourcex=-a0/2;
	removesinkx=bathLength+channelLength+bathLength+a0;
	
	removetopchannely=-b0/2;
	removebottomchannely=channelWidth+b0/2;
	
	
	sourceDensity=(int)(bathLength*bathWidth*sourceBfield/Phi); 
	sinkDensity=(int)(bathLength*bathWidth*sinkBfield/Phi); 
	channelDensity=(int)(channelLength*(channelWidth+b0)*((sourceBfield+sinkBfield)/Phi)/2.0); 	
	
	 
	std::cout << "Tube geometry selected." << std::endl;
	
}

void GeometryTube::InitialiseVortices() const
{

	std::cout << "Initialising Vortices..." << std::endl;
	std::cout << "   " << "sourceDensity: " << sourceDensity << std::endl;
	std::cout << "   " << "sinkDensity: " << sinkDensity << std::endl;
	
	
	std::cout << "   " << pos_file_name << std::endl;
	
	std::ifstream myfile (pos_file_name.c_str());
	
	if (myfile.is_open()) // Get all particle positions from file
	{
		std::cout << "   " << "Initial Vortex Positions From File" << std::endl;
		
		double xval;
		double yval;
		char type;
		
		while ( myfile.good() )
		{
			myfile >> type >> xval >> yval;
						
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			newVortex.set_type(type);
			if (type=='A') AParticles.push_back(newVortex);  
			else if OtherParticles.push_back(newVortex);
				
	
		}
		myfile.close();
	
	}
	else // Make random mobile particles and CE particles
	{
		InitialiseRandomMobileParticles();
		InitialiseCEParticles()
	
	}
	
	std::cout << "   " << "initialiseVortices() created " << vorticesList->size() << " vortices." << std::endl << std::endl;
}
        
void InitialiseRandomAParticles()
{         
		std::cout << "   " << "no start data" << std::endl;
		
		for (int i = 0; i<(sourceDensity);i++)
		{
			
			double xval = bathLength*(rand() % 1000)/1000.0;
			double yval = bathWidth*(rand() % 1000)/1000.0;
	
			CParticle newVortex;
	
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList.push_back(newVortex);
	
		}
		for (int i = 0; i<(sinkDensity);i++)
		{
			CParticle newVortex;
			double xval= bathLength*(rand() % 1000)/1000.0+bathLength+channelLength;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList.push_back(newVortex);
		}	
		for (int i = 0; i<(channelDensity);i++)
		{
			double xval = channelLength*(rand() % 1000)/1000.0+bathLength;
			double yval = channelWidth*(rand() % 1000)/1000.0;
			
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList.push_back(newVortex);
			
		}	         
}         
             
void GeometryTube::ReplaceEscapedVortices() const
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList.begin();
			p!=AParticlesList.end(); ++p)
	{
		
		double x = p->get_x();
		double y = p->get_y();
		
		if (x <= removesourcex || 
				x >= removesinkx )
		{		
			//std::cout << "Escaped particle at (" <<  x << ", " << y << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
		
			if (x <= removesourcex)  // put back in source
			{
				double xval = bathLength*(rand() % 1000)/1000.0;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				p->set_pos(xval,yval);
			}
			else if  (x >= removesinkx) // put back in sink 
			{
				double xval= bathLength*(rand() % 1000)/1000.0+channelLength+bathLength;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				p->set_pos(xval,yval);
			}
		}	
	
	  //continue;
		if ((y <= removetopchannely) &&  (x > removesourcex) )
		{
			//std::cout << "wrapped at (" <<  x << ", " << y << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p->set_pos(x,y+channelWidth+b0);
		}
		else if ((y >= removebottomchannely) &&  (x > removesourcex) )
		{
			//std::cout << "wrapped at (" <<  x << ", " << y << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p->set_pos(x,y-channelWidth-b0);
		}
	  
	
	
	}
	
	
	
	
}

void GeometryTube::InitialiseCEParticles()
{
	std::cout << "Initialising CE..." << std::endl;
	
	double locala0=a0;
	double localb0=b0;
		
	xlo=-3*a0;
	ylo=-10*localb0;
	
	xhi = 2*bathLength+channelLength+3*locala0;
	yhi = channelWidth+10*localb0-0.1*b0+3*localb0/2.0;
	
	double xPos;
	double yPos=ylo;
	
		
	while (yPos<yhi)
	{
		xPos=xl;

	    while (xPos < xhi)
	    {
			CParticle newPin;
			newPin.set_pos(xPos,yPos+localb0/2.0);
			newVortex.set_type('W');
			OtherParticlesList.push_back(newPin);
					
			newPin.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
			newVortex.set_type('W');
			AParticlesList.push_back(newPin);
			
			xPos=xPos+locala0;
			
		}
		
		yPos=yPos+2*localb0;
	}
	
	yhi = yhi + 3*localb0/2.0; //fiddle
		
	//etch source, sink and channel
	bool removed;
	
	std::list<CParticle>::iterator p= pinsList->begin();
	
	while (p!=pinsList->end())
	{
		
		removed=false;
		
		if (  p->get_x() > etchsourcex && p->get_x() < etchsinkx
		  )
		{
			p=pinsList->erase(p);
			removed=true;
		}
		
		if (removed==false) { ++p; }
	
			
	}
	
	std::cout << "   initialisePins() created " << pinsList->size() << " 'CE' vortices." << std::endl <<std::endl;
	
	
	
}

void GeometryTube::AddParticlesForDT(std::list<CParticle> & list_) const
{
	newList.insert(AParticlesList.begin(),AParticlesList.end());
	
	for (std::list<CParticle>::iterator p = AParticlesList.begin();
			p!=AParticlesList.end(); ++p )
	{
		if (p->get_y() <= 2*b0)
		{
			CParticle newVortex;
			newVortex = (*p);
			newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelWidth+b0);
			newList.push_back(newVortex);
		}
		else if (p->get_y() >= channelWidth-2*b0)
		{
			CParticle newVortex;
			newVortex = (*p);
			newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelWidth-b0);
			newList.push_back(newVortex);
		}
	}
}

double GeometryTube::GetRemovalSourceX() const
{
	return 5*a0;
} 

double GeometryTube::GetRemovalSinkX() const
{
	return 2*bathLength+channelLength - 5*a0;
} 

void GeometryTube::UpdateBathDensities() const
{
	int t = sim.get_t();
	
	static double lastchangedSource=0;
	static double lastchangedSink=0;
	
	double actualSource = calcSourceB();
	double actualSink = calcSinkB();
	
	if (t%100==0) std::cout << sourceBfield << "(" << actualSource << ")       (" << sinkBfield << "(" << actualSink << ")" << std::endl;
		
	double expectedSource = sourceBfield;
	double expectedSink = sinkBfield;
	
	
	// calculate zone densities
	int sourceCount=0;
	int sinkCount=0;
	
	// relaxation time
	double relaxation_time = a0/fabs(sim.get_tAvSAvVelX())/channelWidth*b0*.75; 
	
	if (t%1000==0) std::cout << "relaxation time: " << relaxation_time << std::endl;
	
	// count densities
	for (std::list<CParticle>::iterator p = AParticlesList.begin();
			p != AParticlesList.end(); ++p)
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
		}
		if (actualSource>expectedSource)
		{
			if (RemoveParticleFromBath("source")) lastchangedSource=t;
		}
	}
	
	if (dt*t-dt*lastchangedSink>=relaxation_time)		// update sink after relaxation time
	{
		if (actualSink<expectedSink)
		{
			if (AddParticleToBath("sink")) lastchangedSink=t;
		}
		if (actualSink>expectedSink)
		{
			if (RemoveParticleFromBath("sink")) lastchangedSink=t;
		}
	}
	
}

bool GeometryTube::AddParticleToBath(std::string location_) const
{
		int t = sim.get_t();
	
		if (location_.compare("source") != 0 && location_.compare("sink") != 0)
			throw std::runtime_error("AddParticleToBath() location_ is not source or sink");

		if (location_.compare("source") == 0)
		{
		
			CParticle newVortex;
				
			double xval = (2*a0)*(rand() % 1000)/1000.0;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			xval=xval+a0/2.0;
		
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList.push_back(newVortex);
				
		}
		else if (location_.compare("sink") == 0)
		{

			CParticle newVortex;
				
			double xval = (2*a0)*(rand() % 1000)/1000.0;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			xval=xval-a0/2.0;
				
			// offset to end of system
			xval = 2*bathLength + channelLength-xval;
				
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList.push_back(newVortex);
				
		}		
		
		if (location_.compare("source") == 0)
		{
			std::cout << t << "  +                    " << std::endl;	
		}	
		if (location_.compare("sink") == 0)
		{
			std::cout << t << "                      +" << std::endl;	
		}	
		
		return true;
}

bool GeometryTube::RemoveParticleFromBath(std::string location_) const
{
		int t = sim.get_t();
		
		if (location_.compare("source") != 0 && location_.compare("sink") != 0)
			throw std::runtime_error("RemoveParticleFromBath() location_ is not source or sink");
		
		double removalx;
		
		int vortex_to_remove = -1;
		
		// make a vector of pointers to particles in the target bath
			
		std::vector<CParticle> targetVortices;
		std::vector<CParticle> otherVortices;
		
		if (location_.compare("source") == 0)
		{
		
			removalx = GetRemovalSourceX();
		
			// calculate which vortices are in removal zone
			
			for (std::list<CParticle>::iterator p = AParticlesList.begin();
				p != AParticlesList.end(); ++p) {
				
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
			removalx = GetRemovalSinkX();
			
		  // calculate which vortices are in removal zone
			
			for (std::list<CParticle>::iterator p = AParticlesList.begin();
				p != vorticesList->AParticlesList.end(); ++p) {
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
			vorticesList->clear();
			
			std::copy( otherVortices.begin(), otherVortices.end(), std::back_inserter( *AParticlesList ) );
			std::copy( targetVortices.begin(), targetVortices.end(), std::back_inserter( *AParticlesList ) );
			
			if (location_.compare("source") == 0)
			{
				std::cout << t << "  -                    " << std::endl;	
			}	
			if (location_.compare("sink") == 0)
			{
				std::cout << t << "                      -" << std::endl;	
			}	
				
			
			return true;
		}
		
	return false;
		
}

double GeometryTube::calcSinkB() const
{
	
	double aaverage=0;
	int numa=0;
	for (std::list<CDelLine>::iterator p = triangulatedLinesList->begin();
				p!=triangulatedLinesList->end(); ++p)
	{
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if ( (midx > bathLength+channelLength-binsize/2.0 && midx<bathLength+channelLength+binsize/2.0) &&
				(midy>0 && midy<channelWidth))
		{
		
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));		
			aaverage=aaverage+linelength;
			numa++;	
		
		}				
	}
	aaverage=aaverage/(double)numa;
	
	return 2*Phi/(sqrt((double)3)*aaverage*aaverage);	 // B effective
	
}

double GeometryTube::calcSourceB() const
{
	
	double aaverage=0;
	int numa=0;
	for (std::list<CDelLine>::iterator p = triangulatedLinesList->begin();
				p!=triangulatedLinesList->end(); ++p)
	{
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
		
		double ymin=0;
		double ymax=channelWidth;				
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
	
	return 2*Phi/(sqrt((double)3)*aaverage*aaverage);	 // B effective
	
	
}

void GeometryTube::PerStepAnalysis()
{
	  OutputVortexPositions(); 
}

void GeometryTube::EndofSimAnalysis()
{
	OutputFinalVortexPositions();
	OutputPinsList();
	OutputAverages();
	CalculateAndOutputNd();
	CalculateAndOutputAvVel();
	
}

void GeometryTube::PerStepUpdates()
{
	UpdateBathDensities();
	ReplaceEscapedVortices();	
	
}

void GeometryTube::WrapVortices(std::list<CParticle>& vorticesList_) const
{
	
	double wrapsize = forcerange;
	double ysize = channelWidth+b0;
	std::list<CParticle> wrappedVorticesList;
	wrappedVorticesList=vorticesList_;
	 
	for (std::list<CParticle>::iterator p = vorticesList_.begin();
		p!=vorticesList_.end(); ++p )
	{
			if (p->get_y() <= wrapsize)  // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+ysize);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			if (p->get_y() >= ysize-wrapsize) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-ysize);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
	}
	/*std::cout << "VorticesList---------" << std::endl;
	for (std::list<CParticle>::iterator p = vorticesList_.begin();
		p!=vorticesList_.end(); ++p )
	{
		std::cout << p->get_x() << " " << p->get_y() << std::endl;
	}
	std::cout << "---------------" << std::endl;
	*/
	vorticesList_=wrappedVorticesList;
	
	/*std::cout << "wrappedVorticesList---------" << std::endl;
	for (std::list<CParticle>::iterator p = vorticesList_.begin();
		p!=vorticesList_.end(); ++p )
	{
		std::cout << p->get_x() << " " << p->get_y() << std::endl;
	}
	std::cout << "---------------" << std::endl;
	*/
}

void GeometryTube::CalculateAndOutputAvVel()
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

void GeometryTube::OutputFinalVortexPositions()
{
	// At the end of the simulation, output vortex positions
	if (t==simulation_time+1)
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
		
		std::cout << "Writing final vortex positions...done" << std::endl;
	}
}

void GeometryTube::OutputCE()
{
/*	static bool PinsOutputDone = false;
	// Output the pinsList	
	if (PinsOutputDone==true)
			return;
	
	PinsOutputDone=true;
	for (std::list<CParticle>::iterator p = OtherParticlesList.begin();
				p!=OtherParticlesList.end();++p) {
		*fileOutputter.getFS("partcilesfile") << " " << p->get_x() << " " << p->get_y() << std::endl;
	
	}
	
	std::cout << "Writing CE pins positions...done" << std::endl;*/
}

void GeometryTube::OutputParticlePositions()
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

void GeometryTube::OutputAverages()
{	
	if (simulation_time+1!=t) throw std::runtime_error("Averages must be output at the end of the simulation.");
	
		*fileOutputter.getFS("avfile") << "Time and space averaged quantities" << std::endl
					 << "  velx of channel vortices: " << avXVel/t << std::endl
					 << "  vely of channel vortices: " << avYVel/t << std::endl
					 << "  M2 (just stochastic term): " << get_M2Average() << std::endl
					 << "  M2 (all terms): " << get_M2FullAverage() << std::endl;
	
	std::cout << "Writing final averages...done" << std::endl;
	
}

void GeometryTube::CalculateAndOutputNd()
{
	if (t%triangulationInterval!=0)
		return;
		
	Nd=0;
	Nv=0;
	Nmis=0;
	for (std::list<CParticle>::iterator p = delVortexList.begin();
		p!=delVortexList.end(); ++p)
	{
		if (p->get_x() <get_bathLength() || p->get_x()>get_bathLength()+get_channelLength()) continue;
			if (false==p->get_ghost())
			{
				Nv++;
				
				if (6!=p->get_coord_num())
				{
					Nmis++;
				}
				
			} 
			
	}
	
	Nd=Nmis/(double)Nv;
	*fileOutputter.getFS("Nd") << t << " " << Nd << std::endl;
}

bool isA (const CParticle & i) { return i.get_Type()=='A' ? true : false; }

void GeometryTube::GetIParticles(std::list<CParticle> & iList) const
{
	// just A particles
	iList.insert(vorticesList.begin(),vorticesList.end());
	iList.remove_if(isA);
	
} 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
