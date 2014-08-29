#include "CFileOutput.hpp"

#include <sstream>
#include <iterator>
#include <iostream>
#include <memory>
#include <iomanip>

#if defined(__WINDOWS__)
  #include <windows.h>
#else
  #include <sys/stat.h>
#endif


	
CFileOutput::CFileOutput()
{
	jobDirectory="";
	fsptr= new std::ofstream*[MAXFILES];
	
}


//*************************************************************************************************************
// 
//	Sets Job Directory
//
//
//*************************************************************************************************************
	
void CFileOutput::setJobDirectory(std::string jobNum_)
{
	filecount=0;
	jobDirectory=jobNum_;
	
	// create directory		
	#if defined (__WINDOWS__)
		CreateDirectory (jobDirectory.c_str(), NULL);
	#else
		mkdir(jobDirectory.c_str(),0755);	
	#endif
	

}

//*************************************************************************************************************
// 
//	Close all files
//
//
//*************************************************************************************************************
	
CFileOutput::~CFileOutput()
{
		//for (fsMap::iterator p = files.begin();
		//		p!=files.end();++p)
		//{
		//	(*(p->second)).close();
			
		//}
		
		for(int i=0;i<filecount;++i)
			newFile[i].close();

}


//*************************************************************************************************************
// 
//	add file stream
//
//
//*************************************************************************************************************
	
void CFileOutput::addFileStream(std::string streamName_, std::string fileName_)
{
	std::ostringstream oss;
	
	// file name
	oss << jobDirectory << "//" << fileName_;
	
	
	newFile[filecount].open(oss.str().c_str());
	newFile[filecount].precision(5);
	//newFile[filecount].setf( std::ios::scientific, std:: ios::floatfield );	
	
	fsptr[filecount]=&newFile[filecount];
	
	files.insert ( std::pair<std::string, std::ofstream* >(streamName_, fsptr[filecount] ));
	std::cout << "Initialised " << streamName_ << "(" << fileName_ << ")" << std::endl;
	
	
	filecount++;
}


//*************************************************************************************************************
// 
//	getFS
//
//		returns a pointer to the filestream
//
//*************************************************************************************************************
		
std::ofstream * CFileOutput::getFS(std::string streamName_)
{		
		return files.find(streamName_)->second;
	
}

//*************************************************************************************************************
// 
//	isFS
//
// 		Checks if there is a valid FS
//
//*************************************************************************************************************
	
bool CFileOutput::isFS(std::string streamName_)
{		
		return files.find(streamName_)!=files.end();
}
