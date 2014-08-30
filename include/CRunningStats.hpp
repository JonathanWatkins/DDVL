#ifndef CRUNNINGSTATS_HPP
#define CRUNNINGSTATS_HPP

// B. P. Welford algorithm for numerically accurate running mean and variance.
// Presented in Donald Knuth's Art of Computer Programming, Vol 2, page 232, 3rd edition

#include <stdexcept>
#include <cmath>

class CRunningStats
{
	int k;
	
	double Mk;
	double Mkminus1;
	
	double Sk;
	double Skminus1;
	
	double Mk2;
	
public:
	CRunningStats()
	{
		k=0;
	
		Mk=0;
		Mkminus1=0;
		
		Mk2=0;
		
		Sk=0;
		Skminus1=0;
	
	}
	
	~CRunningStats()
	{}
	
	void add(double x_)
	{
		++k;
		
		/*if (1==k)
    {
			Mkminus1 = Mk = x_;
			Skminus1 = 0.0;
    }
    else
		{
			Mk = Mkminus1 + (x_ - Mkminus1)/(double)k;
			Sk = Skminus1 + (x_ - Mkminus1)*(x_ - Mk);

			// set up for next iteration
			Mkminus1 = Mk; 
			Skminus1 = Sk;
    }*/
    
    Mk+=x_;
    Mk2+=x_*x_;
    
		
	}
	
	double get_mean() const
	{
		//return (k > 0) ? Mk : 0.0;
		return (k>0) ? Mk/k : 0.0;
	}
	
	double get_se() const
	{
		//return ( (k > 1) ? Mk/(double)(k - 1) : 0.0 );
		double radix = Mk2 - Mk*Mk/k;
    
    if(radix < 0) throw std::runtime_error("CRunningStats:get_se():  Radix is negative");
    
    return std::sqrt(radix)/k;
	
	}
	
	int get_numDataPoints() const
	{
		return (int)k;
	}

	void reset()
	{
		k=0;
		Mk=0;
		Mk2=0;
	}

};

#endif
