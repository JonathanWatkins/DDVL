#ifndef CDELLINE_HPP
#define CDELLINE_HPP

#include <cmath>

class CDelLine {
    double x1, y1, x2, y2,length;
    
public:
  
  CDelLine(double x1_, double y1_, double x2_, double y2_)
  : x1(x1_), y1(y1_), x2(x2_), y2(y2_) 
  {
		length=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	};
  
	CDelLine() {
		x1=0;
		y1=0;
		x2=0;
		y2=0;
		length=0;
		
	};
	
	~CDelLine() {};
	
	
	 
    void set_points (double valx1, double valy1, double valx2, double valy2)	{
			//if (valx1<valx2) {
				x1=valx1;
				y1=valy1;
				x2=valx2;
				y2=valy2;
				length = std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
			/*}
			else {
				x1=valx2;
				y1=valy2;
				x2=valx1;
				y2=valy1;
			}*/
			
			
	};
	
	bool operator ==(CDelLine const a) const
	{
		return 	(std::fabs(x1-a.x1)<0.00001 && std::fabs(y1 - a.y1) < 0.00001 && std::fabs(x2 - a.x2) < 0.00001 && std::fabs(y2 - a.y2) < 0.00001) ||
				(std::fabs(x2 - a.x1) < 0.00001 && std::fabs(y2 - a.y1) < 0.00001 && std::fabs(x1 - a.x2) < 0.00001 && std::fabs(y1 - a.y2) < 0.00001);
		
	//	return x < a.x;
	};
	
	/*bool operator !=(CCoord const a) const
	{
		return x != a.x || y != a.y;
		
	//	return x < a.x;
	};*/
	
		double get_x1 ()	{
			return x1;
	};
    
    	double get_y1 ()	{
			return y1;
	};
	
		double get_x2 ()	{
			return x2;
	};
	
		double get_y2 ()	{
			return y2;
	}
    double get_length ()	{
			return length;
	}
  
    
};
#endif
