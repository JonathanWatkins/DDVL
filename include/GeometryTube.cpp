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

GeometryTube::GeometryTube(CSimulation * sim_)
{
	sim=new CSimulation;
	sim=sim_;
	
	triangulatedParticlesList = new std::list<CParticle>;
    triangulatedLinesList = new std::list<CDelLine>;
    AParticlesList = new std::list<CParticle>;
    OtherParticlesList = new std::list<CParticle>;
	
	bathLength = 0;
	bathWidth = 0;
	channelLength = 0;
	channelWidth = 0;
	sourceBfield = 0;
	sinkBfield = 0;
	a0 = 0;
	b0 = 0;
	pos_file_name = "";
	Phi = 0;
	forcerange=0;
	dt=0;
	
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
        
    Nd=0;
	Nv=0;
	Nmis=0;    
	
}

GeometryTube::~GeometryTube()
{
	delete triangulatedParticlesList;
    delete triangulatedLinesList;
    delete AParticlesList;
    delete OtherParticlesList;
	
}

void GeometryTube::LoadBatchFile()
{
	std::cout << "Loading job batch file..." << std::endl;
	
	// geometry variables
	sim->ReadVariableFromBatchFile(a0, "GeneralParameters.a0");
	b0=(std::sqrt((double)3)/2.0)*a0;
		
	sim->ReadVariableFromBatchFile(channelLength, "Geometry.channelLength");
	sim->ReadVariableFromBatchFile(channelWidth, "Geometry.channelWidth");
	sim->ReadVariableFromBatchFile(sourceBfield, "Geometry.sourceBfield");
	sim->ReadVariableFromBatchFile(sinkBfield, "Geometry.sinkBfield");
	sim->ReadVariableFromBatchFile(bathLength, "Geometry.bathLength");
	sim->ReadVariableFromBatchFile(bathWidth, "Geometry.bathWidth");
	sim->ReadVariableFromBatchFile(Phi, "Interactions.Phi");
	sim->ReadVariableFromBatchFile(forcerange, "GeneralParameters.forceRange");
	sim->ReadVariableFromBatchFile(dt, "GeneralParameters.dt");
	
	channelLength*=a0;
	channelWidth*=b0;
	bathLength*=a0;
	channelWidth*=b0;
	
	
	// analysis variables
	sim->ReadVariableFromBatchFile(binsize, "GeneralParameters.binSize");
	bool posfile;
	sim->ReadVariableFromBatchFile(posfile, "InputData.altPosFile");
	if (posfile == true)
	{
			sim->ReadVariableFromBatchFile(pos_file_name, "InputData.altPosFileName");
	
	}
	
				
			
	// interactions
	
	
	
	
	
	std::cout << "   Job Header loaded.\n\n";
	
	
}

void GeometryTube::InitialiseGeometry()
{
	LoadBatchFile();
	InitialiseParameters();
	InitialiseVortices();
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
	

	std::cout << "a0: " << a0
			  << "Phi: " << Phi
			  << "channelLength: " << channelLength
			  << "channelWidth: " << channelWidth
			  << "bathLength: " << bathLength
			  << "bathWidth: " << bathWidth
			  << "sourceBfield: " << sourceBfield
			  << "sinkBfield: " << sinkBfield
			  << std::endl;
			   
	 
	std::cout << "Tube geometry selected." << std::endl;
	
}

void GeometryTube::InitialiseVortices()
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
			if (type=='A') AParticlesList->push_back(newVortex);  
			else OtherParticlesList->push_back(newVortex);
				
	
		}
		myfile.close();
	
	}
	else // Make random mobile particles and CE particles
	{	
		InitialiseRandomMobileParticles();
		InitialiseCEParticles();
	
	}
	
	std::cout 	<< "   " << "initialiseVortices() created " << AParticlesList->size() << " Langevin vortices" 
				<< " and " << OtherParticlesList->size() << " other votices." << std::endl << std::endl;
}
        
void GeometryTube::InitialiseRandomMobileParticles()
{         
		std::cout << "   " << "no start data" << std::endl;
		
		for (int i = 0; i<(sourceDensity);i++)
		{
			
			double xval = bathLength*(rand() % 1000)/1000.0;
			double yval = bathWidth*(rand() % 1000)/1000.0;
	
			CParticle newVortex;
	
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList->push_back(newVortex);
	
		}
		for (int i = 0; i<(sinkDensity);i++)
		{
			CParticle newVortex;
			double xval= bathLength*(rand() % 1000)/1000.0+bathLength+channelLength;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList->push_back(newVortex);
		}	
		for (int i = 0; i<(channelDensity);i++)
		{
			double xval = channelLength*(rand() % 1000)/1000.0+bathLength;
			double yval = channelWidth*(rand() % 1000)/1000.0;
			
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			newVortex.set_type('A');
			AParticlesList->push_back(newVortex);
			
		}	         
}         
             
void GeometryTube::ReplaceEscapedVortices()
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p!=AParticlesList->end(); ++p)
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
		xPos=xlo;

	    while (xPos < xhi)
	    {
			CParticle newVortex;
			newVortex.set_pos(xPos,yPos+localb0/2.0);
			newVortex.set_type('W');
			OtherParticlesList->push_back(newVortex);
					
			newVortex.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
			newVortex.set_type('W');
			OtherParticlesList->push_back(newVortex);
			
			xPos=xPos+locala0;
			
		}
		
		yPos=yPos+2*localb0;
	}
	
	yhi = yhi + 3*localb0/2.0; //fiddle
		
	//etch source, sink and channel
	bool removed;
	
	std::list<CParticle>::iterator p= OtherParticlesList->begin();
	
	while (p!=OtherParticlesList->end())
	{
		
		removed=false;
		
		if (  p->get_x() > etchsourcex && p->get_x() < etchsinkx
		  )
		{
			p=OtherParticlesList->erase(p);
			removed=true;
		}
		
		if (removed==false) { ++p; }
	
			
	}
	
}

void GeometryTube::AddParticlesForDT(std::list<CParticle> & iList)
{
	// Add A particles
	iList.clear();
	std::list<CParticle>::iterator it = iList.end();
	iList.insert(it,AParticlesList->begin(),AParticlesList->end());
	
	WrapVortices(iList);
	
	static bool once = false;
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
	
}

double GeometryTube::GetRemovalSourceX() const
{
	return 5*a0;
} 

double GeometryTube::GetRemovalSinkX() const
{
	return 2*bathLength+channelLength - 5*a0;
} 

void GeometryTube::UpdateBathDensities()
{
	int t = sim->get_t();
	
	static double lastchangedSource=0;
	static double lastchangedSink=0;
	
	double actualSource = CalcSourceB();
	double actualSink = CalcSinkB();
	
	if (t%sim->get_framedataInterval()==0) std::cout << sourceBfield << "(" << actualSource << ")       (" << sinkBfield << "(" << actualSink << ")" << std::endl;
		
	double expectedSource = sourceBfield;
	double expectedSink = sinkBfield;
	
	
	// calculate zone densities
	int sourceCount=0;
	int sinkCount=0;
	
	// relaxation time
	double relaxation_time = a0/fabs(avXVel/t)/channelWidth*b0*.75; 
	
	if (t%(sim->get_framedataInterval()*10)==0) std::cout << "relaxation time: " << relaxation_time << std::endl;
	
	// count densities
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
			p != AParticlesList->end(); ++p)
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

bool GeometryTube::AddParticleToBath(std::string location_)
{
		int t = sim->get_t();
	
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
			AParticlesList->push_back(newVortex);
				
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
			AParticlesList->push_back(newVortex);
				
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

bool GeometryTube::RemoveParticleFromBath(std::string location_)
{
		int t = sim->get_t();
		
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
			
			for (std::list<CParticle>::iterator p = AParticlesList->begin();
				p != AParticlesList->end(); ++p) {
				
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
			
			for (std::list<CParticle>::iterator p = AParticlesList->begin();
				p != AParticlesList->end(); ++p) {
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
			AParticlesList->clear();
			
			std::list<CParticle>::iterator it = AParticlesList->end();
			AParticlesList->insert(it, otherVortices.begin(), otherVortices.end());
			AParticlesList->insert(it, targetVortices.begin(), targetVortices.end());
			
			//std::copy( otherVortices.begin(), otherVortices.end(), std::back_inserter( AParticlesList ) );
			//std::copy( targetVortices.begin(), targetVortices.end(), std::back_inserter( AParticlesList ) );
			
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

double GeometryTube::CalcSinkB() const
{
	
	/*double aaverage=0;
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
	*/
	
	int np=0;
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
				p!=AParticlesList->end(); ++p)
	{
		double x = p->get_x();
		double y = p->get_y();
						
		if ( (x > bathLength+channelLength-binsize/2.0 && x<bathLength+channelLength+binsize/2.0) &&
				(y>-b0/2 && y<channelWidth+b0/2))
		{
		
			np++;
		}				
	}
	
	double B= Phi*np/binsize/channelWidth;
	
	//std::cout << "Sink B: "  << B  << " " << 2*Phi/(sqrt((double)3)*aaverage*aaverage) << std::endl;
	
	return B;
	
}

double GeometryTube::CalcSourceB() const
{
	
	/*double aaverage=0;
	int numa=0;
	
	//if (triangulatedLinesList->size() == 0) throw std::runtime_error("GeometryTube::CalcSourceB() No triangulated particles");
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
	*/
	
	int np=0;
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
				p!=AParticlesList->end(); ++p)
	{
		double x = p->get_x();
		double y = p->get_y();
						
		if ( (x > bathLength-binsize/2.0 && x<bathLength+binsize/2.0) &&
				(y>-b0/2 && y<channelWidth+b0/2))
		{
		
			np++;
		}				
	}
	
	double B= Phi*np/binsize/channelWidth;
	
	//std::cout << "Sink B: "  << B  << " " << 2*Phi/(sqrt((double)3)*aaverage*aaverage) << std::endl;
	
	return B;
	
	
}

void GeometryTube::PerStepAnalysis()
{
	  OutputParticlePositions(); 
	  CalculateAndOutputAvVel();
}

void GeometryTube::EndofSimAnalysis()
{
	OutputFinalParticlePositions();
	OutputAverages();
	CalculateAndOutputNd();
	
	
}

void GeometryTube::PerStepUpdates()
{
	UpdateBathDensities();
	ReplaceEscapedVortices();	
	
}

void GeometryTube::WrapVortices(std::list<CParticle>& iList)
{
		
	// Add periodic y particles	
	double wrapsize = forcerange;
	double ysize = channelWidth+b0;
	
	for (std::list<CParticle>::iterator p = AParticlesList->begin();
		p!=AParticlesList->end(); ++p )
	{
			if (p->get_y() <= wrapsize)  // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+ysize);
				newVortex.set_type('B');
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
			if (p->get_y() >= ysize-wrapsize) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-ysize);
				newVortex.set_type('B');
				newVortex.set_ghost();
				iList.push_back(newVortex);
			}
	}

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
	for(std::list<CParticle>::iterator p = AParticlesList->begin();
	    p != AParticlesList->end(); p++)
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
	int t = sim->get_t();
	*sim->get_FS("framevel") << t << " " << spaceSumX/double(count) << " " << spaceSumY/double(count) << " " << avXVel/t << " " << avYVel/t << std::endl;
	if (t%sim->get_framedataInterval()==0) std::cout << t << " " << spaceSumX/double(count) << " " << spaceSumY/double(count) << " " << avXVel/t << " " << avYVel/t << std::endl;
	
	
}

void GeometryTube::OutputFinalParticlePositions()
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

void GeometryTube::OutputParticlePositions()
{
	
	int t = sim->get_t();	
	if (t==1) *sim->get_FS("guifile") << "# This file contains frame data\n"
																	<< "# { t, numofparticles, {x1,y1,velx1,vely1,coordnum1},...,{xN,yN,velxN,velyN,coordnumN}}" << std::endl; 
	
	if (t%sim->get_triangulationInterval()!=0 || t%sim->get_framedataInterval()!=0) return;
	
	// counts number of active particles
	
	*sim->get_FS("guifile") << "{" << t << ", " << triangulatedParticlesList->size() << ", ";
	
	
	bool first=true;
	for (std::list<CParticle>::iterator p = triangulatedParticlesList->begin();
			p != triangulatedParticlesList->end(); ++p)
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

void GeometryTube::OutputAverages()
{	
	int t = sim->get_t();
	if (sim->get_simulation_time()+1!=t) throw std::runtime_error("Averages must be output at the end of the simulation.");
	
		*sim->get_FS("avfile") << "Time and space averaged quantities" << std::endl;
		*sim->get_FS("avfile") << "  velx of channel vortices: " << avXVel/t << std::endl;
		*sim->get_FS("avfile") << "  vely of channel vortices: " << avYVel/t << std::endl;
		*sim->get_FS("avfile") << "  M2 (just stochastic term): " << sim->get_M2Average() << std::endl;
		*sim->get_FS("avfile") << "  M2 (all terms): " << sim->get_M2FullAverage() << std::endl;
	
	std::cout << "Writing final averages...done" << std::endl;
	
}

void GeometryTube::CalculateAndOutputNd()
{
	if (sim->get_t()%sim->get_triangulationInterval()!=0)
		return;
		
	Nd=0;
	Nv=0;
	Nmis=0;
	for (std::list<CParticle>::iterator p = triangulatedParticlesList->begin();
		p!=triangulatedParticlesList->end(); ++p)
	{
		if (p->get_x() <bathLength || p->get_x()>bathLength+channelLength) continue;
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
	*sim->get_FS("Nd") << sim->get_t() << " " << Nd << std::endl;
}

std::list<CParticle> * GeometryTube::GetIParticles()
{
	return AParticlesList;
}

void GeometryTube::GetJParticles(std::list<CParticle> & iList)
{
	// Add A particles
	iList.clear();
	std::list<CParticle>::iterator it = iList.end();
	iList.insert(it,AParticlesList->begin(),AParticlesList->end());


	// Add CE particles
	iList.insert(it,OtherParticlesList->begin(),OtherParticlesList->end());

	WrapVortices(iList);

}
 
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
