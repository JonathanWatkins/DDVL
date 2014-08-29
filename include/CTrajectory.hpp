#ifndef CTRAJECTORY_HPP
#define CTRAJECTORY_HPP

#include <vector>

class CTrajectory
{

private:
	std::vector<double> xlist;
	std::vector<double> ylist;
	int NumPoints;

    
    
public:
	CTrajectory()
	{
		NumPoints=0;
	}
	
	~CTrajectory() {}
	
	void add_TrajectoryPoint(const double x_, const double y_)
	{
				xlist.push_back(x_);
				ylist.push_back(y_);
				NumPoints++;
		
	}
	
	std::vector<double>* get_xlist()
	{
		return &xlist;
	}
	
	std::vector<double>* get_ylist()
	{
		return &ylist;
		
	}
	
	
	int get_NumPoints() const
	{
		return NumPoints;
	}	  
    
};


#endif
