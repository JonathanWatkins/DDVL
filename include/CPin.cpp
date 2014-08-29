#include "CPin.hpp"

CPin::CPin()
{
	x=0;
	y=0;
	A=0;
	r=0;
	
}
	 
CPin::CPin(double x_, double y_, double A_, double r_) 
: x(x_), y(y_), A(A_), r(r_)
{}
 
double CPin::get_x() const
{
	return x;
}

double CPin::get_y() const
{
	return y;
}

double CPin::get_A() const
{
	return A;
}

double CPin::get_r() const
{
	return r;
}

void CPin::set_pos(double x_, double y_)
{
	x=x_;
	y=y_;
}

void CPin::set_A(double A_)
{
	A=A_;
}

void CPin::set_r(double r_)
{
	r=r_;
}

void CPin::set_x(double x_)
{
	x=x_;
}

void CPin::set_y(double y_)
{
	y=y_;
}
