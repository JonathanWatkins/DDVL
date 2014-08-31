//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryChannel.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "GeometryChannel.hpp"
#include "CSimulation.hpp"
#include "CParticle.hpp"


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryChannel::GeometryChannel(CSimulation & sim_)
:		sim(sim_)
,   vorticesList(sim_.get_vorticesList())
,		pinsList(sim_.get_pinsList())
,   delLinesList(sim_.get_delLinesList())
,		Phi(sim_.get_Phi())
,		a0(sim_.get_a0())
,		b0(sim_.get_b0())
,		dt(sim_.get_dt())
,		binsize(sim.get_binsize())
,   pos_file_name(sim_.GetPosFileName())
,		pins_file_name(sim_.GetPinsFileName())
,		jobBatchFileLocation(sim_.get_jobBatchFileLocation())
{

	LoadBatchFile();
	
	
	// calculate system parameters
	channelOffset = (bathWidth-channelWidth)/2.0;
	
	etchsourcex0 =0; 
	etchsourcey0= -channelOffset;
	etchsourcex1 =bathLength-0.1*a0; 
	etchsourcey1= channelWidth+channelOffset;
	
	etchsinkx0 =bathLength+channelLength+0.1*a0; 
	etchsinky0= -channelOffset;
	etchsinkx1 =bathLength+channelLength+bathLength; 
	etchsinky1= channelWidth+channelOffset;
	
	etchchannelx0=0;
	etchchannely0=0;
	etchchannelx1=bathLength+channelLength+bathLength;
	etchchannely1=channelWidth;
	
	removesourcex=-a0/2;
	removesinkx=bathLength+channelLength+bathLength+a0;
	removechannelx0=bathLength+a0/2;
	removechannelx1=bathLength+channelLength-a0/2;
	removesourcey0=removetopchannely=channelOffset-3.0*b0/2;
	removesourcey1=removebottomchannely=channelWidth+channelOffset+3.0*b0/2;
	
	sourceDensity=(int)(bathLength*bathWidth*sourceBfield/Phi); 
	sinkDensity=(int)(bathLength*bathWidth*sinkBfield/Phi); 
	channelDensity=(int)(channelLength*(channelWidth+b0)*((sourceBfield+sinkBfield)/Phi)/2.0); 	
	 
	
}

void GeometryChannel::LoadBatchFile()
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
	else {
	
		std::cout << "   " << "no start data" << std::endl;
		//file = false;
	
		for (int i = 0; i<(sourceDensity);i++) {
				
			double xval,yval;
			
			xval = bathLength*(rand() % 1000)/1000.0;
			yval = bathWidth*(rand() % 1000)/1000.0;
			//yval = yval-channelOffset;
			
			CParticle newVortex;
				
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
		
		}
		for (int i = 0; i<(sinkDensity);i++)
		{
		
			CParticle newVortex;
			double xval= bathLength*(rand() % 1000)/1000.0+channelLength+bathLength;
			double yval = bathWidth*(rand() % 1000)/1000.0;
			//yval = yval-channelOffset;
			//cout << "sink vortex: " << xval << ", " << yval << endl;
			newVortex.set_pos(xval,yval);
			vorticesList->push_back(newVortex);
		}	

		for (int i = 0; i<(channelDensity);i++)
		{
			
			double xval = channelLength*(rand() % 1000)/1000.0+bathLength;
			double yval = channelWidth*(rand() % 1000)/1000.0;
			//cout << "channel vortex: " << xval << ", " << yval << endl;
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			
			vorticesList->push_back(newVortex);
		
		}	
		
	}
   
	std::cout << "   " << "initialiseVortices() created " << vorticesList->size() << " vortices." << std::endl << std::endl;
	
}
         
void GeometryChannel::ReplaceEscapedVortices() const
{
	
 	for (std::list<CParticle>::iterator p = vorticesList->begin();
			p!=vorticesList->end(); ++p)
	{
 	
		if (p->get_x() <= removesourcex || p->get_y() <= removesourcey0 || 
				p->get_x() >= removesinkx  || p->get_y() >= removesourcey1 )
		{		
			std::cout << "Escaped particle at (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
		
			if (p->get_x() < bathLength)  // put back in source
			{
				double xval = bathLength*(rand() % 1000)/1000.0;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				p->set_pos(xval,yval);
			}
			else if  (p->get_x() > bathLength + channelLength) // put back in sink 
			{
				double xval= bathLength*(rand() % 1000)/1000.0+channelLength+bathLength;
				double yval = bathWidth*(rand() % 1000)/1000.0;
				p->set_pos(xval,yval);
			}
			else  // put back in channel
			{
				double xval = channelLength*(rand() % 1000)/1000.0+bathLength;
				double yval = channelWidth*(rand() % 1000)/1000.0;
				p->set_pos(xval,yval);
			}
		}	
	}
	
	
}

void GeometryChannel::InitialisePins()
{
	
	//double offset=channelOffset;//channelWidth/2.0;
	
	std::cout << "Initialising pins..." << std::endl;
	
	
	double locala0=a0;
	double localb0=b0;
	double xPos;
	
	
	firstPin.set_pos(-3*a0,-6*localb0-channelOffset);
		
	std::cout << "   " << pins_file_name << std::endl;

	std::ifstream myfile (pins_file_name.c_str());

		
	if (myfile.is_open()) 
	{
		std::cout << "   " << "Initial CE Positions From File" << std::endl;
		
		double xval;
		double yval;
	
		while ( myfile.good() )
		{
			myfile >> xval;
			myfile >> yval;
		
			CParticle newVortex;
			newVortex.set_pos(xval,yval);
			pinsList->push_back(newVortex);
	
		}
		myfile.close();
	
	}
	else
	{
		
		double yPos=firstPin.get_y();
		
		//std::cout << "Bottom Limit: " << yPos << std::endl; 
		//std::cout << "Top Limit: " << bathWidth+6*localb0-channelOffset << std::endl; 
		
		while (yPos<bathWidth+6*localb0-channelOffset-0.1*b0)
		{
			xPos=firstPin.get_x();
		
			while (xPos<bathLength+channelLength+bathLength+3*locala0)
			{
		
				CParticle newPin;
				newPin.set_pos(xPos,yPos+localb0/2.0);
				pinsList->push_back(newPin);
		
				newPin.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
				pinsList->push_back(newPin);
				
				xPos=xPos+locala0;
		
			}
		
			yPos=yPos+2*localb0;
			std::cout << "   yPos: " << yPos << std::endl;
		}
		
		//etch source, sink && channel
		bool removed;
		std::list<CParticle>::iterator p= pinsList->begin();
		
		while (p!=pinsList->end())
		{
			removed=false;
			
				//etch source
				if (  p->get_x() > etchsourcex0 && p->get_x() < etchsourcex1
					&& p->get_y() > etchsourcey0 && p->get_y() < etchsourcey1 )
				{
					p=pinsList->erase(p);
					removed=true;
			
				}
				else if (  p->get_x() > etchchannelx0 && p->get_x() < etchchannelx1
					&& p->get_y() > etchchannely0 && p->get_y() < etchchannely1 )
				{
					p=pinsList->erase(p);
					removed=true;
				
				}
				else if (  p->get_x() > etchsinkx0 && p->get_x() < etchsinkx1
					&& p->get_y() > etchsinky0 && p->get_y() < etchsinky1 )
				{
					p=pinsList->erase(p);
					removed=true;
				}
			
		
			if (removed==false) { ++p; }
		
		
		}
	}
	std::cout << "   initialisePins() created " << pinsList->size() << " CE vortices." << std::endl <<std::endl;
	
	
	
}

void GeometryChannel::AddParticlesForDT(std::list<CParticle> & vorticesList_) const
{

	for(std::list<CParticle>::iterator p=pinsList->begin();
				p!=pinsList->end();p++)
		{
			if (  p->get_x() > 0 && p->get_x() < 2*bathLength+channelLength
					&& p->get_y() > etchchannely0-a0 && p->get_y() < etchchannely1+a0 )
			{
				(*p).set_ghost();	
				vorticesList_.push_back(*p);
			}
		}	
	
	
}

void GeometryChannel::WrapSystem() const
{}

void GeometryChannel::InitialiseDisorder() const
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

double GeometryChannel::GetRemovalSourceX() const
{
	return 5*a0;
} 

double GeometryChannel::GetRemovalSinkX() const
{
	return 2*bathLength+channelLength - 5*a0;
} 

CParticle GeometryChannel::GetFirstPin() const
{
	return firstPin;
}



void GeometryChannel::UpdateBathDensities() const
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

bool GeometryChannel::AddParticleToBath(std::string location_) const
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

bool GeometryChannel::RemoveParticleFromBath(std::string location_) const
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


double GeometryChannel::calcSinkB() const
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

double GeometryChannel::calcSourceB() const
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

void GeometryChannel::WrapVortices(std::list<CParticle>& vorticesList_) const
{
	return;
	/* no wrap
	double wrapsize = sim.get_forceRange();
	double ysize = channelWidth;
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
			else if (p->get_y() >= ysize-wrapsize) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-ysize);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
	}
	
	vorticesList_=wrappedVorticesList;
	*/
}


 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
