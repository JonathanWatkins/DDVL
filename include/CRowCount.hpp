#ifndef CROWCOUNT_HPP
#define CROWCOUNT_HPP


class CRowCount {
    double x;
		int count;
    
public:
  
  CRowCount(double x_, int count_)
  : x(x_), count(count_)
  {};
  
	CRowCount() {
		x=0;
		count=0;
	};
	
	~CRowCount() {};
	
	void set_val (double x_, int count_)
  {
		x=x_;
		count=count_;
	};
	
	double get_x() const
	{
		return x;
	};
	
	int get_count() const
	{
		return count;
	};
	
	
    
};
#endif
