//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	GeometryChannel.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <stdexcept>
#include <list>
#include <iterator>
#include <iostream>

#include "GeometryChannel.h"
#include "CSimulation.hpp"
#include "CParticle.hpp"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

GeometryChannel::GeometryChannel(CSimulation & sim_)
:   vorticesList(sim_.get_vorticesList())
,   removesourcex(sim_.get_removesourcex())
,   removesourcey0(sim_.get_removesourcey0())
,   removesinkx(sim_.get_removesinkx())
,   removesourcey1(sim_.get_removesourcey1())
,   removechannelx0(sim_.get_removechannelx0())
,   removechannelx1(sim_.get_removechannelx1())
,   removetopchannely(sim_.get_removetopchannely())
,   removebottomchannely(sim_.get_removebottomchannely())
,   t(&sim_.get_t())
{}
                

void GeometryChannel::RemoveEscapedVortices() const
{

	std::list<CParticle>::iterator p = vorticesList.begin();
 
 	while (p != vorticesList.end())
	{
		bool removed=false;
		
		if (p->get_x() <= removesourcex)
		{
			
			//OutputTrajectory(p);
			
			std::cout << "Removed(type1) at " << *t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		}
		else if (p->get_y() <= removesourcey0)
		{
			
			//OutputTrajectory(p);
			
			std::cout << "Removed(type2) at " << *t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
		}
		else if (p->get_x() >= removesinkx)
		{
			
			//OutputTrajectory(p);
			
			std::cout << "Removed at(type3) " << *t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} 
		else if (p->get_y() >= removesourcey1)
		{
			
			//OutputTrajectory(p);
			
			std::cout << "Removed at(type4) " << *t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
		} 
		else if ( (p->get_x() >= removechannelx0 && p->get_x()<= removechannelx1) 
				&& ( p->get_y() <= removetopchannely  || p->get_y() >= removebottomchannely ) )
		{
			//OutputTrajectory(p);
			
			std::cout << "Channel Removed at " << *t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << std::endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} 
		
		if (removed==false) { ++p; }
	
	}
	
	
}

void GeometryChannel::InitialisePins() const
{}

void GeometryChannel::InitialiseVortices() const
{}

void GeometryChannel::AddParticlesForDT() const
{}

void GeometryChannel::WrapSystem() const
{}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
