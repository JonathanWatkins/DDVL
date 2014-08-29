#include "Utilities.hpp"

#include <boost/math/special_functions.hpp>

namespace Utilities
{

	const double pi = acos(-1.0);

}

void Utilities::wrapVortices(std::list<CParticle>& vorticesList_, double ysize_, double wrapsize_)
{
	
	std::list<CParticle> wrappedVorticesList;
	wrappedVorticesList=vorticesList_;
	
	for (std::list<CParticle>::iterator p = vorticesList_.begin();
		p!=vorticesList_.end(); ++p )
	{
			if (p->get_y() <= wrapsize_)  // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= ysize_-wrapsize_) //channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
	}
	
	vorticesList_=wrappedVorticesList;
	
}
	

void Utilities::wrapVorticesPeriodic(std::list<CParticle>& vorticesList_, double xsize_, double ysize_, double wrapsize_)
{
	std::list<CParticle> wrappedVorticesList;
	wrappedVorticesList=vorticesList_;
	
	for (std::list<CParticle>::iterator p = vorticesList_.begin();
		p!=vorticesList_.end(); ++p )
	{
		// wrap vortices on tube
		//if (p->get_x() >=0) {
			if (p->get_y() <= wrapsize_) // forcerange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+ysize_); // channelWidth
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= ysize_- wrapsize_) // channelWidth-forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			if (p->get_x() <= wrapsize_)  // forceRange
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+xsize_,newVortex.get_y()); // channelLength
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_x() >= xsize_ - wrapsize_)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-xsize_,newVortex.get_y()); // channelLengtg
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			
			// corners
			if (p->get_y() <= wrapsize_ && p->get_x() <=wrapsize_)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+xsize_,newVortex.get_y()+ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= ysize_-wrapsize_ && p->get_x() >= xsize_-wrapsize_)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-xsize_,newVortex.get_y()-ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() <= wrapsize_ && p->get_x() >= xsize_-wrapsize_)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()-xsize_,newVortex.get_y()+ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= ysize_-wrapsize_ && p->get_x() <= wrapsize_)
			{
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x()+xsize_,newVortex.get_y()-ysize_);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			
			
			
			
			
			
			
		//}
		
		/*else if (p->get_x() < 0) {  // wrap vortices on cone
			if (p->get_y() <= -fabs(p->get_x())*tan(funnelAngle)+forceRange) {
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelHeight+b0+2.0*fabs(p->get_x())*tan(funnelAngle) );
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= channelHeight+fabs(p->get_x())*tan(funnelAngle)-forceRange) {
				CParticle newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelHeight-b0-2.0*fabs(p->get_x())*tan(funnelAngle));
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			
		}*/
	}
	
	vorticesList_=wrappedVorticesList;

}


double Utilities::gaussianRand()
{
	double result;
	do
	{
		double u1 = (float)rand()/(float)RAND_MAX;
		double u2 = (float)rand()/(float)RAND_MAX;
		result = std::sqrt((double)-2.0*std::log(u1))*std::cos(2.0*pi*u2);
	}
	while ( boost::math::isnan((double)result)==true || boost::math::isinf(result)==true  || result >1 || result <-1);
	
	if (result != result)
	{
		std::cout << "ERR: gaussian rand result: " << "(" << result << ")" << std::endl;
	}
	
	if (boost::math::isinf(result)) std::cout << "inf: gaussian rand result: " << "(" << result << ")" << std::endl;
		 
	return result;
}
