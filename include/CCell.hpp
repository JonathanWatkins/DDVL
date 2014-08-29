#ifndef CCELL_HPP
#define CCELL_HPP

#include <vector>
#include <list>
#include "CParticle.hpp"
#include <iostream>
#include <string>
#include <iterator>

class CCell
{
	public:
	
	std::list<CParticle> cellList;

	public:
		
	CCell(){};
	
	~CCell(){};
	
	std::list<CParticle> * get_cellList();
	
	void add_particle(CParticle vortex_);
	
	void list_particles();
	
	void set_cellList(std::list<CParticle> cellList_);
	
	CCell(const CCell& lhs_)
	{
		
		//std::cout << "Copy-construct called" << std::endl;
		std::copy( (lhs_.cellList).begin(), (lhs_.cellList).end(), std::back_inserter( cellList ) );
		
	};
	
	CCell& operator= (const CCell& lhs_)
	{
		
		//std::cout << "Copy-assign called" << std::endl;
		std::copy( (lhs_.cellList).begin(), (lhs_.cellList).end(), std::back_inserter( cellList ) );
		
		return *this;
	};
	
	void clearlist() { cellList.clear(); };
	
};


#endif
