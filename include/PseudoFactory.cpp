//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	PseudoFactory.cpp
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <stdexcept>
#include "PseudoFactory.h"

#include "Input.h"

#include "OptionBase.h"
#include "EuroPut.h"
#include "EuroCall.h"
#include "OptionAverageRate.h"

#include "ProcessBase.h"
#include "GBM_process.h"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	interface
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

PseudoFactory::PseudoFactory()
{
    inp_ = new Input;
}
  
PseudoFactory::~PseudoFactory()
{
    delete inp_;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	getter
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double PseudoFactory::GetS0() const {return inp_->GetS0();}
double PseudoFactory::Getr() const {return inp_->Getr();}
double PseudoFactory::Getsig() const {return inp_->Getsig();}
char PseudoFactory::GetPtype() const {return inp_->GetPtype();}

double PseudoFactory::GetX() const {return inp_->GetX();}
double PseudoFactory::GetT() const {return inp_->GetT();}
char PseudoFactory::GetOtype() const {return inp_->GetOtype();}
        
long PseudoFactory::GetM() const {return inp_->GetM();}
long PseudoFactory::GetN() const {return inp_->GetN();}
  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	CreateOption
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

OptionBase * PseudoFactory::CreateOption()
{
    char o_type = inp_->GetOtype();
    
    switch(o_type)
    {
        case 'c':  return new EuroCall(*this);            break;
        case 'p':  return new EuroPut(*this);             break;
        case 'a':  return new OptionAverageRate(*this);   break;
        default:   throw std::runtime_error("Valuation::CreateOption:  Bad character");
    }
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	CreateProcess
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ProcessBase * PseudoFactory::CreateProcess()
{
    char p_type = inp_->GetPtype();
    
    switch(p_type)
    {
        case 'g':  return new GBM_process(*this);         break;
        default:   throw std::runtime_error("Valuation::CreateProcess:  Bad character");
    }
}
 
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//	end
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
