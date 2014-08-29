#ifndef CPOSITIONBASE_HPP
#define CPOSITIONBASE_HPP

class CPositionBase
{
	public:
		
	 //CPositionBase(){};
	 
	 //CPositionBase(double x_, double y_, double A_, double r_) ;
	 
	virtual ~CPositionBase(){}
	 
	virtual double get_x() const=0;

	virtual double get_y() const=0;

	virtual void set_pos(double , double )=0;
	
	virtual void set_x(double x_) =0;
	
	virtual void set_y(double y_)=0;
	
};


#endif
