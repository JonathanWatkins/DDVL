#ifndef CVERSION_HPP
#define CVERSION_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

class CVersion
{
	private:
	int major;
	int minor;
	int maintenance;
	
	public:
		
	CVersion()
	{
		major=0;
		minor=0;
		maintenance=0;
	};
	
	CVersion(int major_, int minor_, int maintenance_)
	{};
	
	~CVersion(){};
	
	void set_versionStr(std::string version_)
	{
		
		std::ostringstream oss;
		
		
		major= (int)atoi(version_.substr(0,1).c_str());
		
		
		minor= (int)atoi(version_.substr(02,1).c_str());
		
		
		maintenance= (int)atoi(version_.substr(4,1).c_str());
	
	}
	
	std::string get_versionStr() const
	{
		std::ostringstream oss;
		
		oss << major << "." << minor << "." << maintenance;
		
		return oss.str();
	}
		
	bool operator< (const CVersion& rhs_)
	{
		//std::cout << get_versionStr() << " < " << rhs_.get_versionStr() << std::endl;
		int thisver=100*major+10*minor+maintenance;
		int rhsver =100*rhs_.major+10*rhs_.minor+rhs_.maintenance;
		if (thisver<rhsver)
			return true;
			
		return false;
	};
	
	bool operator<= (const CVersion& rhs_)
	{
		int thisver=100*major+10*minor+maintenance;
		int rhsver =100*rhs_.major+10*rhs_.minor+rhs_.maintenance;
		if (thisver<=rhsver)
			return true;
			
		return false;
	};
	
	
	
	bool operator== (const CVersion& rhs_)
	{
		return (rhs_.major == major && rhs_.minor == minor && rhs_.maintenance==maintenance);
	
	};
	
	
};


#endif
