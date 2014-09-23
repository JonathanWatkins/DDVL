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
:		sim(sim_)
,   vorticesList(sim_.get_vorticesList())
,		pinsList(sim_.get_pinsList())
,   delLinesList(sim_.get_delLinesList())
,		Phi(sim_.get_Phi())
,		a0(sim_.get_a0())
,		b0(sim_.get_b0())
,		dt(sim_.get_dt())
,		binsize(sim.get_binsize())
,   forcerange(sim.get_forceRange())
,   pos_file_name(sim_.GetPosFileName())
,		pins_file_name(sim_.GetPinsFileName())
,		jobBatchFileLocation(sim_.get_jobBatchFileLocation())
{

	LoadBatchFile();
	
	// calculate system parameters
	channelOffset = (bathWidth-channelWidth)/2.0;
	
	etchsourcex = 0; 
	etchsinkx = bathLength+channelLength+bathLength; 
	
	removesourcex=-a0/2;
	removesinkx=bathLength+channelLength+bathLength+a0;
	
	removetopchannely=-b0/2;//channelOffset-3.0*b0/2;
	removebottomchannely=channelWidth+b0/2;//channelWidth+channelOffset+3.0*b0/2;
	
	
	sourceDensity=(int)(bathLength*bathWidth*sourceBfield/Phi); 
	sinkDensity=(int)(bathLength*bathWidth*sinkBfield/Phi); 
	channelDensity=(int)(channelLength*(channelWidth+b0)*((sourceBfield+sinkBfield)/Phi)/2.0); 	
	 
	std::cout << "Tube geometry selected." << std::endl;
}

void GeometryTube::LoadBatchFile()
{
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(jobBatchFileLocation, pt);
	
	channelLength=pt.get<double>("Geometry.channelLength")*a0;
	channelWidth=pt.get<double>("Geometry.channelWidth")*b0;
	sourceBfield=pt.get<double>("Geometry.sourceBfield");
	sinkBfield=pt.get<double>("Geometry.sinkBfield");
	bathLength=pt.get<double>("Geometry.bathLength")*a0;
	bathWidth=pt.get<double>("Geometry.bathWidth")*b0;
	
}

void GeometryTube::InitialiseVortices() const
{

	std::cout << "Initialising Vortices..." << std::endl;
	std::cout << "   " << "sourceDensity: " << sourceDensity << std::endl;
	std::cout << "   " << "sinkDensity: " << sinkDensity << std::endl;
	
	
	std::cout << "   " << pos_file_name << std::endl;
	
	std::ifstream myfile (pos_file_name.c_str());
	
	
	if (myfile.is_open()) 
	{
		std::cout << "   " << "Initial Vortex Positions From File" << std::endl;
		
		double xval;
		double yval;
	
		while ( myfile.good() )
		{
			myfile >> xval;
			myfile >> yval;
		
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
	
		}
		myfile.close();
	
	}
	else
	{
		std::cout << "   " << "no start data" << std::endl;
		
		for (int i = 0; i<(sourceDensity);i++)
		{
			
			double xval = bathLength*(rand() % 1000)/1000.0;
			double yval = bathWidth*(rand() % 1000)/1000.0;
	
			CParticle newVortex;
	
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
	
		}
		for (int i = 0; i<(sinkDensity);i++)
		{
			CParticle newVortex;
			double xval= bathLength*(rand() % 1000)/1000.0+bathLength+channelLength;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
		
		
		}	
	
	
	
		for (int i = 0; i<(channelDensity);i++)
		{
			double xval = channelLength*(rand() % 1000)/1000.0+bathLength;
			double yval = channelWidth*(rand() % 1000)/1000.0;
			
			
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
			
		}	
	
	}
	
	std::cout << "   " << "initialiseVortices() created " << vorticesList->size() << " vortices." << std::endl << std::endl;
}
         
void GeometryTube::ReplaceEscapedVortices() const
{
	// replaces particles that escape the source and wraps particles in y direction along the channel
 	for (std::list<CParticle>::iterator p = vorticesList->begin();
			p!=vorticesList->end(); ++p)
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

void GeometryTube::InitialisePins()
{
	
	std::cout << "Initialising pins..." << std::endl;
	
	double locala0=a0;
	double localb0=b0;
	
	double xPos;
	
		
	firstPin.set_pos(-3*a0,-10*localb0-channelOffset);
	
	double yPos=firstPin.get_y();
	
	//std::cout << "Bottom Limit: " << yPos << std::endl; 
	//std::cout << "Top Limit: " << bathWidth+6*localb0-channelOffset << std::endl; 
	
	while (yPos<channelWidth+10*localb0-channelOffset-0.1*b0)
	{
		xPos=firstPin.get_x();

	    while (xPos<2*bathLength+channelLength+3*locala0)
	    {
			CParticle newPin;
			newPin.set_pos(xPos,yPos+localb0/2.0);
			pinsList->push_back(newPin);
					
			newPin.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
			pinsList->push_back(newPin);
			
			xPos=xPos+locala0;
			
		}
		
		yPos=yPos+2*localb0;
	}
		
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

void GeometryTube::AddParticlesForDT(std::list<CParticle> & vorticesList_) const
{

	for (std::list<CParticle>::iterator p = vorticesList->begin();
			p!=vorticesList->end(); ++p )
	{
		if (p->get_y() <= 2*b0)
		{
			CParticle newVortex;
			newVortex = (*p);
			newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelWidth+b0);
			newVortex.set_ghost();
			vorticesList_.push_back(newVortex);
		}
		else if (p->get_y() >= channelWidth-2*b0)
		{
			CParticle newVortex;
			newVortex = (*p);
			newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelWidth-b0);
			newVortex.set_ghost();
			vorticesList_.push_back(newVortex);
		}
	}
}

void GeometryTube::WrapSystem() const
{}

void GeometryTube::InitialiseDisorder() const
{
	/*std::cout << "Initialising channel disorder..." << std::endl;
	
	std::cout << "   Number of Gaussian pins per a0^2: " << disorderDensity << std::endl;
	
	double numberChannelPins=channelLength*channelWidth*disorderDensity/a0/a0;	
	std::cout << "   Number of Gaussin pins in this channel: " << numberChannelPins << std::endl;
	
	std::cout << "   Disorder strength: " << disorderStrength << std::endl;
	std::cout << "   Disorder range: " << disorderRange << std::endl;
	
	
	for (int i=0; i< numberChannelPins;i++)
	{
		double x=bathLength+channelLength*(rand()/(double)RAND_MAX);
		double y=channelWidth*(rand()/(double)RAND_MAX);
		CParticle newPin;
		newPin.set_pos(x,y);
		disorderList.push_back(newPin);
	}
	
	std::cout << "   Channel disorder initialised." << std::endl << std::endl;
	*/ 
}

double GeometryTube::GetRemovalSourceX() const
{
	return 5*a0;
} 

double GeometryTube::GetRemovalSinkX() const
{
	return 2*bathLength+channelLength - 5*a0;
} 

CParticle GeometryTube::GetFirstPin() const
{
	return firstPin;
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
	for (std::list<CParticle>::iterator p = vorticesList->begin();
			p != vorticesList->end(); ++p)
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
			vorticesList->push_back(newVortex);
				
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
			vorticesList->push_back(newVortex);
				
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
			
			for (std::list<CParticle>::iterator p = vorticesList->begin();
				p != vorticesList->end(); ++p) {
				
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
			
			for (std::list<CParticle>::iterator p = vorticesList->begin();
				p != vorticesList->end(); ++p) {
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
			
			std::copy( otherVortices.begin(), otherVortices.end(), std::back_inserter( *vorticesList ) );
			std::copy( targetVortices.begin(), targetVortices.end(), std::back_inserter( *vorticesList ) );
			
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
	for (std::list<CDelLine>::iterator p = delLinesList->begin();
				p!=delLinesList->end(); ++p)
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
	for (std::list<CDelLine>::iterator p = delLinesList->begin();
				p!=delLinesList->end(); ++p)
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


 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
