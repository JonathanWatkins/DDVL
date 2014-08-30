//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryChannel.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>
#include <fstream>

#include "GeometryChannel.hpp"
#include "CSimulation.hpp"
#include "CParticle.hpp"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryChannel::GeometryChannel(CSimulation & sim_)
:   vorticesList(sim_.get_vorticesList())
,		pinsList(sim_.get_pinsList())
,		bathLength(sim_.get_bathLength())
,		bathWidth(sim_.get_bathWidth())
,		channelLength(sim_.get_channelLength())
,		channelWidth(sim_.get_channelWidth())
,		sourceBfield(sim_.get_sourceBfield())
,		sinkBfield(sim_.get_sinkBfield())
,		Phi(sim_.get_Phi())
,		a0(sim_.get_a0())
,		b0(sim_.get_b0())
,   pos_file_name(sim_.GetPosFileName())
,		pins_file_name(sim_.GetPinsFileName())
{
	
		
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

void GeometryChannel::AddParticlesForDT() const
{

	for(std::list<CParticle>::iterator p=pinsList->begin();
				p!=pinsList->end();p++)
		{
			if (  p->get_x() > 0 && p->get_x() < 2*bathLength+channelLength
					&& p->get_y() > etchchannely0-a0 && p->get_y() < etchchannely1+a0 )
			{
				(*p).set_ghost();	
				vorticesList->push_back(*p);
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

 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
