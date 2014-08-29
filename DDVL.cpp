#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include "CSimulation.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

int main(int argc, char *argv[])
{		
   	try
	{
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
				std::string	jobBatchFileLocation = argv[1];
				sim.Initialise(jobBatchFileLocation);
			}
		}
			
			
		sim.Run();	
			
		
	}
	catch (const boost::property_tree::file_parser_error &e)
    {   
		std::cout 	<< "A boost runtime error was caught by main\n"
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

