#include "CCell.hpp"


std::list<CParticle> * CCell::get_cellList()
{
	return &cellList;
}

void CCell::add_particle(CParticle vortex_)
{
	cellList.push_back(vortex_);
	
}

void CCell::list_particles()
{
	std::cout << "This list has " << cellList.size() << " particles." <<std::endl;
	for (std::list<CParticle>::iterator p=cellList.begin();
			p!=cellList.end(); ++p)
	{
		
		std::cout << p->get_x() << ", " << p->get_y() << std::endl;
		
		
	}
	
}

void CCell::set_cellList(std::list<CParticle> cellList_)
{
	cellList=cellList_;
	
}



