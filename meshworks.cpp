#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation


#include "CSimulation.hpp"
#include "CSDLGLscreen.hpp"
#include <iostream>
#include <SDL.h>
#include "CEvents.hpp"
#include <sstream>
#include <string>
#include <stdexcept>
#include <boost/exception/all.hpp>

#include <boost/throw_exception.hpp>

typedef boost::error_info<struct my_tag,std::string> my_tag_error_info;

int main(int argc, char *argv[])
{
		
		/*
		 * The number of command line parameters depends on the runtype and geometry
		 * 
		 * runtype = 0  requires argc=5
		 *  eg.   meshworks 0 <starting timestep> <0=no triangulation, 1=triangulate> <jobnumber>  
		 * 
		 * runtype = 1 and 2  geometry = 0 (channel) requires argc=12
		 *  eg.   meshworks <1=write 2=write with gui> <output 0=none 1=all> 
		 * 										0 <sourceBfield (T)> <sinkBfield (T)> 
		 * 										<bathLength (a0)> <bathWidth (b0)> <channelLength (a0)>	<channelWidth (b0)> 
		 * 										<numTimeSteps> <tau>  
		 * 
		 * runtype = 1 and 2 geometry = 1 (tube) requires argc=11
		 *  eg.   meshworks <1=write 2=write with gui> <output 0=none 1=all> 
		 * 										1 <sourceBfield (T)> <sinkBfield (T)> 
		 * 										<bathLength (a0)> <channelLength (a0)>	<channelWidth (b0)> 
		 * 										<numTimeSteps> <tau>  * 
		 */
		  
		int runtype; 
		int outputtype;
		int geometry;
		
		double sourceBfield;
		double sinkBfield;
		
		double bathLength;
		double bathWidth;
		double channelLength;
		double channelWidth;
		
		int simulation_time;
		double temp;
		int starting_time_step;
		
		double triangulateReadRun;
		std::string jobnum;
		
		std::string jobBatchFileLocation;
		
   	try
		{
   	CEvents events;
		
    CSimulation sim;
		
		// initialise from jobHeader file
		if(2==argc)
		{
			std::stringstream oss;
			oss << argv[1];
			if ("help"==oss.str())
			{
				std::cout << 	"* The number of command line parameters depends on the runtype and geometry\n"
									<<	"*\n"
									<<  "* Simulation should be run using meshworks job.bat\n"
									<<  "* where job.bat should be replaced with the name of the job file\n\n" 
									<< 	"* To read previously generated job data runtype = 0  requires argc=5\n"
									<<	"*  eg.   meshwork 0 <starting timestep> <0=no triangulation, 1=triangulate> <jobnumber>\n";
									return 1;
		 	}
		 	else
		 	{ 
				jobBatchFileLocation=argv[1];
				sim.initialise(jobBatchFileLocation);
				runtype=sim.get_runtype(); 
			}
		}
		else if(5 == argc && 0==(int)atof(argv[1]))
		{
			runtype= (int)atof(argv[1]);
			starting_time_step= (int)atof(argv[2]);
			if (1==(int)atof(argv[3])) triangulateReadRun=true;
			else triangulateReadRun=false;
			
			jobnum = argv[4];  
			
			// initialise write job
			//if (1==runtype || 2==runtype) sim.initialiseWrite(runtype,outputtype,geometry,sourceBfield,sinkBfield,bathLength,bathWidth,channelLength,channelWidth,simulation_time,temp);
			
			// initialise read job
			std::stringstream oss;
			oss << jobnum;
			
			if (0==runtype)	sim.initialiseRead(runtype,starting_time_step,triangulateReadRun,oss.str());
		}	
		
		CSDLGLscreen window;
		
		if(0==runtype || 2==runtype)
		{	
			window.initialiseSDLandGL(&sim);
		}
		
		do 
		{
			events.doEvents(&sim);
			
			
			sim.do_step();
			
			if (true==sim.get_draw())
			{
				window.drawSystem();
				
			}
		
			sim.next_t();
		
			//std::cout << std::endl;
		} while(sim.is_running());
    
    sim.clean_up();
		}
		catch (const boost::property_tree::file_parser_error &e)
    {   
			std::cout << "A boost runtime error was caught by main\n"
								<< e.what() << std::endl
                << e.message() << std::endl
                << e.line() << std::endl;
    }	
		catch (const std::runtime_error & e)
		{
			std::cout << "A runtime error was caught by main\n" << e.what() << '\n';
		}
		catch (...)
		{
			std::cout << "An unknown error has occured!\n";
		}
		
		
		

    return 0;
    
}

