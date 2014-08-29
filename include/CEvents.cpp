#include "CEvents.hpp"
#include "SDL.h"
#include <iostream>
#include "CSimulation.hpp"

void CEvents::doEvents(CSimulation* sim_) 
{

	while (SDL_PollEvent(&event))
	{
		// Check for the quit message
		if (event.type == SDL_QUIT)
		{
			// Quit the program
			sim_->end_simulation();
			std::cerr << "SDL_QUIT event"<< std::endl;
			return;
		}
		
		if (event.type == SDL_KEYDOWN)
		{
			if (event.key.keysym.sym == SDLK_ESCAPE)
			{
				// Quit the program
				sim_->end_simulation();
				std::cerr << "SDL_QUIT event"<< std::endl;
				return;
			}
			if (event.key.keysym.sym == SDLK_SPACE)
			{
				if (sim_->is_paused()) sim_->unpause_simulation();
				else if (!sim_->is_paused()) sim_->pause_simulation();
			}
			if (event.key.keysym.sym == SDLK_q)
			{
				sim_->zoom_in();
			}
			if (event.key.keysym.sym == SDLK_a)
			{
				sim_->zoom_out();
			}
			if (event.key.keysym.sym == SDLK_o)
			{
				sim_->pan_left();
			}
			if (event.key.keysym.sym == SDLK_p)
			{
				sim_->pan_right();
			}
			
			
			
		}
	}
}
