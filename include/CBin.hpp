#ifndef CBIN_HPP
#define CBIN_HPP

class CBin {
  
private:
	double valSum,val;
	int count;
    
    
public:
	CBin() {
		valSum=0;
		count=0;
		val=0;
	};
	
	~CBin() {};
	
	
	void set_val (double val_)	{
			valSum += val_;
			count++;
			val = valSum/(double)count;
	};
	
	 
    double get_val ()	{
			return val;
	};
   
  int get_count () {
		return count;
	}; 
    
     
    
};

#endif
