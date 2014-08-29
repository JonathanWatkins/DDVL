#ifndef CPIN_HPP
#define CPIN_HPP

class CPin
{
	private:
		double x,y;
		double A;
		double r;

	public:
		
	 CPin();
	 
	 CPin(double x_, double y_, double A_, double r_);
	 
	 ~CPin(){};
	 
	double get_x() const;

	double get_y() const;

	double get_A() const;

	double get_r() const;
	
	void set_pos(double x_, double y_);
	
	void set_A(double A_);
	
	void set_r(double r_);
	
	void set_x(double x_);
	
	void set_y(double y_);
	
};


#endif
