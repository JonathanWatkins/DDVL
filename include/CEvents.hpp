#ifndef CEVENTS_HPP
#define CEVENTS_HPP

#include "SDL.h"
#include "CSimulation.hpp"

class CEvents
{
private:

	SDL_Event event;


public:

	CEvents(){};
	~CEvents(){};

	void doEvents(CSimulation* sim_);


};
#endif
