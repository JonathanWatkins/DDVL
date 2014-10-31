//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// FileOutput.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "FileOutput.hpp"

#include <sstream>
#include <iterator>
#include <iostream>
#include <memory>
#include <iomanip>
#include <stdexcept>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Handles compilation on windows/linux
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#if defined(__WINDOWS__)
  #include <windows.h>
#else
  #include <sys/stat.h>
#endif

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// FileOutput()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
FileOutput::FileOutput()
{
	filecount=0;
		
	jobDirectory="";
	
	newFile = new std::ofstream[MAXFILES];
	
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// SetJobDirectory()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
void FileOutput::SetJobDirectory(const std::string &  jobNum_)
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

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// ~FileOutput()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
FileOutput::~FileOutput()
{
		
		for(int i=0;i<filecount;++i)
		{
			std::cout << "Closing file " <<  i << std::endl;
			newFile[i].close();
		}
		delete [] newFile;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// AddFileStream()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void FileOutput::AddFileStream(const std::string &  streamName_, const std::string &  fileName_)
{
	if (CheckIsFileStream(streamName_))
	{
		std::stringstream oss;
		oss << "FileOutput::AddFileStream() Attempting to add duplicate filestream  (" << streamName_ << ")"; 
		throw std::runtime_error(oss.str());
	}

	std::ostringstream oss;
	
	// file name
	oss << jobDirectory << "//" << fileName_;

	newFile[filecount].open(oss.str().c_str());
	newFile[filecount].precision(5);
	
	
	files.insert ( std::pair<const std::string , int >(streamName_, filecount));
	
	std::cout << "Initialised " << streamName_ << "(" << fileName_ << ")" << std::endl;
	
	filecount++;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// GetFS()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX	

int FileOutput::GetFS(const std::string &  streamName_)
{		
		fsMap::iterator p = files.find(streamName_);
		
		if (!CheckIsFileStream(streamName_))
		{
			std::stringstream oss;
			oss << "FileOutput::getFS() No such file was created (" << streamName_ << ")"; 
			throw std::runtime_error(oss.str());
		}
		
		return p->second;
	
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// IsFS()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

bool FileOutput::IsFS(const std::string &  streamName_)
{		
		return files.find(streamName_)!=files.end();
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// RegisterOutput()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

void FileOutput::RegisterOutput(const std::string & streamName_, const std::string & text)
{
	int fileid=GetFS(streamName_);
	
	newFile[fileid] << text << std::endl;
	
	
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// CheckIsFilestream()
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

bool FileOutput::CheckIsFileStream(const std::string & streamName_)
{
	fsMap::iterator p = files.find(streamName_);
		
	if (p==files.end()) return false;
	
	return true;

}
