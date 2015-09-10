//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	BinnedAccumulator.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#pragma warning ( disable : 2586  )  // supresses warning a bug due to icc and boost compilation

#include <iostream>
#include <cmath>

#include "BinnedAccumulator.hpp"


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	constructor
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

BinnedAccumulator::BinnedAccumulator(double xmin, double xmax, double binsize) :
xmin_(xmin),
xmax_(xmax),
binsize_(binsize)
{
	
	numbins_ = (xmax_-xmin_)/binsize_;
	values_.resize(numbins_);// = new std::vector<bin>(numbins_);
        
}

BinnedAccumulator::~BinnedAccumulator()
{
    //delete values_;
}

void BinnedAccumulator::AddValue(double x, double f)
{
	int binnum = std::floor((x-xmin_)/binsize_);
	if (binnum < 0 || binnum > numbins_-1) return;
	
	values_[binnum].xsum+=x;
	values_[binnum].fsum+=f;
	values_[binnum].fsquaredsum+=f*f;
	++values_[binnum].N;
	 	
}

void BinnedAccumulator::PrintBinnedAverages() const
{
	std::cout << "binnum" << " " << "bincentre" << " " << "mean" << " " << "se" << " " << "Nvals" << std::endl; 
	
	for (int i = 0; i < numbins_; ++i)
	{
		double acc_xvals = values_[i].xsum;
		double acc_vals = values_[i].fsum;
		double acc_squs = values_[i].fsquaredsum;
		double N = values_[i].N;
		double mean = acc_vals/N;
		double se = std::sqrt(acc_squs - acc_vals*acc_vals/N)/N;
		double bincentre = xmin_+binsize_/2+i*binsize_; //= acc_xvals/N;
		
		std::cout << i << " " << bincentre << " " << mean << " " << se << " " << N << std::endl; 
		
	}
	
}

void BinnedAccumulator::GetBinnedAverages(std::stringstream & oss)
{
	oss << "binnum" << " " << "bincentre" << " " << "mean" << " " << "se" << " " << "Nvals" << std::endl; 
	
	for (int i = 0; i < numbins_; ++i)
	{
		double acc_xvals = values_[i].xsum;
		double acc_vals = values_[i].fsum;
		double acc_squs = values_[i].fsquaredsum;
		double N = values_[i].N;
		double mean = acc_vals/N;
		double se = std::sqrt(acc_squs - acc_vals*acc_vals/N)/N;
		double bincentre = xmin_+binsize_/2+i*binsize_; //= acc_xvals/N;
		
		oss << i << " " << bincentre << " " << mean << " " << se << " " << N << std::endl; 
		
	}
	
}

void BinnedAccumulator::ClearValues()
{
	for (int i = 0; i < numbins_; ++i)
	{
		values_[i].xsum = 0;
		values_[i].fsum = 0;
		values_[i].fsquaredsum = 0;
		values_[i].N = 0;
		
	}
	
	
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
