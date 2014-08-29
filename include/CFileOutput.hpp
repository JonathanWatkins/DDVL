#ifndef CFILEOUTPUT_HPP
#define CFILEOUTPUT_HPP

#include <map>
#include <string>
#include <fstream>
#include <memory>

#define MAXFILES 100

typedef std::map<std::string, std::ofstream*> fsMap;
//typedef std::map<std::string, int> fsMap;

class CFileOutput
{

	public:
		
		CFileOutput();
	 
		~CFileOutput();
		
		void setJobDirectory(std::string jobNum_);
		
		void addFileStream(std::string streamName_, std::string fileName_);
		
		std::ofstream* getFS(std::string streamName_);
		
		bool isFS(std::string streamName_);
	
	private:
	
		int filecount;
		
		fsMap	files;
		
		std::ofstream newFile[MAXFILES];
		
		std::ofstream **fsptr;
				
		std::string jobDirectory;
		 
};


#endif
