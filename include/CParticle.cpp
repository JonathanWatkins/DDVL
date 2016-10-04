#include "CParticle.hpp"
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <stdexcept>


CParticle::CParticle(double x_, double y_, double velx_, double vely_, int coord_num_)
{
	CParticle();
	x=x_;
	y=y_;
	velx=velx_;
	vely=vely_;
	lastx=x;
	lasty=y;
	coord_num=coord_num_;
		
}

CParticle::CParticle(double x_, double y_, double z_, double velx_, double vely_, double velz_, int coord_num_)	
{
	
	CParticle();
	x=x_;
	y=y_;
	z=z_;
	velx=velx_;
	vely=vely_;
	velz=velz_;
	lastx=x;
	lasty=y;
	lastz=z;
	coord_num=coord_num_;
	
}

CParticle::CParticle()
{
	x=0;
	y=0;
	z=0;
	lastx=0;
	lasty=0;
	lastz=0;
	velx=0;
	vely=0;
	velz=0;
	coord_num=0;
	ghost=false;
	a=0;
	id=rand()*rand()+time(0);
	Jxx=0;
	Jyy=0;
	Jxy=0;
	Jyx=0;
	forcex=0;
	forcey=0;
	force_dx=0;
	force_dy=0;
	force_dz=0;
	force_tx=0;
	force_ty=0;
	force_tz=0;
	TrajectoryNumber=-1;
	Tracked=false;
	in_bubble=false;
	burgers_circuit_center=false;
	N=0;
	velx_sum=0;
	vely_sum=0;
	velz_sum=0;
	//numTrajPoints=0;
	//trajectoryx = new std::vector<double>;
	//trajectoryy = new std::vector<double>;
	type='\0';
	
}

double CParticle::get_x() const
{
	return x;
}
	
double CParticle::get_y() const
{
	return y;
}

double CParticle::get_z() const
{
	return z;
}

double CParticle::get_lastx() const
{
	return lastx;
}

double CParticle::get_lasty() const
{
	return lasty;
}

double CParticle::get_lastz() const
{
	return lastz;
}

double CParticle::get_velx() const
{
	return velx;
}
	
double CParticle::get_vely() const
{
	return vely;
}

double CParticle::get_velz() const
{
	return velz;
}

double CParticle::get_midx() const
{
	return 0.5*(lastx+x);
}

double CParticle::get_midy() const
{
	return 0.5*(lasty+y);
}

double CParticle::get_midz() const
{
	return 0.5*(lastz+z);
}
	
int CParticle::get_coord_num() const
{
	return coord_num;
}
	
void CParticle::set_x(double x_)
{
	x=x_;
}
	
void CParticle::set_y(double y_)
{
	y=y_;
}

void CParticle::set_z(double z_)
{
	z=z_;
}
			
void CParticle::set_pos(double x_, double y_)
{
	x=x_;
	y=y_;
}

void CParticle::set_pos(double x_, double y_, double z_)
{
	x=x_;
	y=y_;
	z=z_;
}

void CParticle::set_lastpos(double lastx_, double lasty_)
{
	lastx=lastx_;
	lasty=lasty_;

}

void CParticle::set_lastpos(double lastx_, double lasty_, double lastz_)
{
	lastx=lastx_;
	lasty=lasty_;
	lastz=lastz_;

}	
	
/*void CParticle::set_velx(double velx_)
{
	velx=velx_;
	velx_sum+=velx;
    N++;
}
	
void CParticle::set_vely(double vely_)
{
	vely=vely_;
	vely_sum+=vely;

}*/
	
void CParticle::set_vel(double velx_, double vely_)
{
	velx=velx_;
	vely=vely_;
	velx_sum+=velx;
	vely_sum+=vely;
    N++;	
}

void CParticle::set_vel(double velx_, double vely_, double velz_)
{
	velx=velx_;
	vely=vely_;
	velz=velz_;
	velx_sum+=velx;
	vely_sum+=vely;
	velz_sum+=velz;
    N++;
	
}

	
void CParticle::set_coord_num(int coord_num_)
{
	coord_num=coord_num_;
}

double CParticle::get_a () const
{
	return a;
}
    
void CParticle::set_a (double a_) 
{
	a=a_;
}

bool CParticle::get_ghost() const 
{
	return ghost;
}
		
void CParticle::set_ghost() 
{
	ghost = true;
}

void CParticle::set_Jyy(double Jyy_)
{
	Jyy=Jyy_;
}

double CParticle::get_Jyy() const
{
	return Jyy;
}

void CParticle::set_Jxx(double Jxx_) 
{
	Jxx=Jxx_;
}

double CParticle::get_Jxx() const
{
	return Jxx;
}

void CParticle::set_Jxy(double Jxy_) 
{
	Jxy=Jxy_;
}

double CParticle::get_Jxy() const
{
	return Jxy;
}

void CParticle::set_Jyx(double Jyx_)
{
	Jyx=Jyx_;
}

double CParticle::get_Jyx() const
{
	return Jyx;
}

int CParticle::get_id() const
{
	return id;
}

double CParticle::get_forcex() const
{
	return forcex;
}

double CParticle::get_forcey() const
{
	return forcey;
}

double CParticle::get_forcez() const
{
	return forcez;
}

void CParticle::set_force(double forcex_, double forcey_)
{
	forcex=forcex_;
	forcey=forcey_;
	
}

void CParticle::set_force(double forcex_, double forcey_, double forcez_)
{
	forcex=forcex_;
	forcey=forcey_;
	forcez=forcez_;
	
}

void CParticle::set_in_bubble()
{
	in_bubble=true;
}

bool CParticle::get_in_bubble() const
{
	return in_bubble;
}

void CParticle::coordPlusOne()
{
	coord_num++;
}

double CParticle::get_velx_mean()
{
    return N==0 ? 0 : velx_sum/N;

}

double CParticle::get_vely_mean()
{
    return N==0 ? 0 : vely_sum/N;
	
}

double CParticle::get_velz_mean()
{
    return N==0 ? 0 : velz_sum/N;
	
}

/*void CParticle::add_trajectoryPoint()
{
	(*trajectoryx).push_back(x);
	(*trajectoryy).push_back(y);
	numTrajPoints++;
	
}

std::vector<double> * CParticle::get_trajectoryx() const
{
		return trajectoryx;
}

std::vector<double> * CParticle::get_trajectoryy() const
{
		return trajectoryy;
}

int CParticle::get_numTrajPoints() const
{
	return numTrajPoints;
}*/

int CParticle::get_TrajectoryNumber() const
{
	if (TrajectoryNumber==-1) throw std::runtime_error("CParticle::get_TrajectoryNumber() Trajectory Number should not return -1");
	
	return TrajectoryNumber;
	
}

void CParticle::set_TrajectoryNumber( const int & TrajectoryNumber_)
{
		TrajectoryNumber=TrajectoryNumber_;
}

void CParticle::set_Tracked()
{
	Tracked=true;
}

bool CParticle::get_Tracked() const
{
	return Tracked;
}
