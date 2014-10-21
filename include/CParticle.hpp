#ifndef CPARTICLE_HPP
#define CPARTICLE_HPP

#include  <vector>
#include <iostream>

class CParticle 
{
	
private:
	
	double x,y,lastx, lasty, velx,vely,a,forcex,forcey;
	int coord_num;
	bool ghost, Tracked, in_bubble, burgers_circuit_center;
	int id;
  double Jxy,Jyx,Jxx,Jyy;
  
  double force_dx, force_dy;
  double force_tx, force_ty;
  
  char type;
  
   
  //std::vector<double> *trajectoryx;
  //std::vector<double> *trajectoryy;
  
  //int numTrajPoints;
  int TrajectoryNumber;  // the element of the trajectories Vector
  
   
   
    
public:
		
	CParticle();
	
	CParticle(double x_, double y_, double velx_, double vely_, int coord_num_);
	
	/*CParticle(double x_, double y_, 
		double lastx_, double lasty_,
		double velx_, double vely_, 
		double coord_num_,
		double forcex_, double forcey_,
		double ghost_,
		int id_,
		double Jxy_, double Jyx_, double Jxx_, double Jyy_)
		:
		x(x_), y(y_), 
		lastx(lastx_), lasty(lasty_),
		velx(velx_), vely(vely_), 
		coord_num(coord_num_),
		forcex(forcex_), forcey(forcey_),
		ghost(ghost_),
		id(id_),
		Jxy(Jxy_), Jyx(Jyx_), Jxx(Jxx_), Jyy(Jyy_)
		
		{};
	*/
	
	void copyAssign_NoTraj(const CParticle & lhs_)
	{
		x = lhs_.x;
		y = lhs_.y;
		lastx = lhs_.lastx;
		lasty = lhs_.lasty;
		velx = lhs_.velx;
		vely = lhs_.vely;
		coord_num = lhs_.coord_num;
		a = lhs_.a;
		forcex = lhs_.forcex;
		forcey = lhs_.forcey;
		ghost = lhs_.ghost;
		id = lhs_.id;
		Jxy = lhs_.Jxy;
		Jyx = lhs_.Jyx;
		Jxx = lhs_.Jxx;
		Jyy = lhs_.Jyy;
		
		force_dx=lhs_.force_dx;
		force_dy=lhs_.force_dy;
		force_tx=lhs_.force_tx;
		force_ty=lhs_.force_ty;
		
		TrajectoryNumber= lhs_.TrajectoryNumber;
		Tracked=lhs_.Tracked;
		in_bubble=lhs_.in_bubble;
		burgers_circuit_center=lhs_.burgers_circuit_center;
		type=lhs_.type;
		
		
	}
	
	
	
	
	~CParticle()
	{
		//	delete trajectoryx;
		//	delete trajectoryy;
			
			
	};
	
	double get_x() const;
	
	double get_y() const;
	
	double get_midx() const;
	
	double get_midy() const;
	
	double get_velx() const;
	
	double get_vely() const;
	
	int get_coord_num() const;
	
	void set_x(double x_);
	
	void set_y(double y_);
		
	void set_pos(double x_, double y_);
	
	void set_lastpos(double lastx_, double lasty_);
	
	void set_velx(double velx_);
	
	void set_vely(double vely_);
	
	void set_vel(double velx_, double vely_);
	
	void set_coord_num(int coord_num_);
	
	double get_lastx() const;
	
	double get_lasty() const;
		
	double get_a () const;

	void set_a (double a_); 

	bool get_ghost() const;
		
	void set_ghost(); 

	void set_Jyy(double Jyy_);

	double get_Jyy() const;

	void set_Jxx(double Jxx_); 
	
	double get_Jxx() const;
	
	void set_Jxy(double Jxy_);
	
	double get_Jxy() const;
	
	void set_Jyx(double Jyx_);
	
	double get_Jyx() const;
	
	int get_id() const;
	
	double get_forcex() const;
	
	double get_forcey() const;
	
	void set_force(double forcex_, double forcey_);
	
	void coordPlusOne();
	
	bool get_in_bubble() const;
	
	void set_in_bubble();
	
	int get_TrajectoryNumber() const; 													// return the element number of the particlesTrajectory vector
	
	void set_TrajectoryNumber( const int & TrajectoryNumber_);	// sets the trajectory number
	
	void set_burgers_circuit_center() { burgers_circuit_center=true;}
	
	bool get_burgers_circuit_center() const { return burgers_circuit_center;}
		
	void set_type(char type_) { type = type_; }
	
	char get_type() const { return type; }

		
	//void add_trajectoryPoint();
	
	//std::vector<double> * get_trajectoryx() const;
	
	//std::vector<double> * get_trajectoryy() const;
	
	//int get_numTrajPoints() const;
	
	
	bool operator< (const CParticle& lhs_)
	{
		if(lhs_.get_x() == x)
				return lhs_.get_y() >= y;
		else if(lhs_.get_x() != x)
				return lhs_.get_x() >= x;
			
		return false;
	};
	
	bool operator== (const CParticle& lhs_)
	{
		return (lhs_.get_x() == x && lhs_.get_y() == y);
				
	};
	
/*	CParticle(const CParticle& lhs_)
	{
		//static int timesCalled=0;
		//timesCalled++;
		//std::cout << "Copy-constructor called " << timesCalled << std::endl;
		
		x = lhs_.x;
		y = lhs_.y;
		lastx = lhs_.lastx;
		lasty = lhs_.lasty;
		velx = lhs_.velx;
		vely = lhs_.vely;
		coord_num = lhs_.coord_num;
		a = lhs_.a;
		forcex = lhs_.forcex;
		forcey = lhs_.forcey;
		ghost = lhs_.ghost;
		id = lhs_.id;
		Jxy = lhs_.Jxy;
		Jyx = lhs_.Jyx;
		Jxx = lhs_.Jxx;
		Jyy = lhs_.Jyy;
		
		//trajectoryx = new std::vector<double>;
		//trajectoryy = new std::vector<double>;
		
		//numTrajPoints=lhs_.numTrajPoints;
		//if (numTrajPoints>1) std::cout << "copied particle with " << numTrajPoints << " trajectories." << std::endl;
		
		for (std::vector<double>::iterator p=(*(lhs_.trajectoryx)).begin();
				p!=(*(lhs_.trajectoryx)).end(); ++p)
		{
			(*trajectoryx).push_back(*p);
		}
		
		for (std::vector<double>::iterator p=(*(lhs_.trajectoryy)).begin();
				p!=(*(lhs_.trajectoryy)).end(); ++p)
		{
			(*trajectoryy).push_back(*p);
		}
		
		
				
	};
	
	CParticle& operator= (const CParticle& lhs_)
	{
		
		std::cout << "Copy-assign called" << std::endl;
		
		
		x = lhs_.x;
		y = lhs_.y;
		lastx = lhs_.lastx;
		lasty = lhs_.lasty;
		velx = lhs_.velx;
		vely = lhs_.vely;
		coord_num = lhs_.coord_num;
		a = lhs_.a;
		forcex = lhs_.forcex;
		forcey = lhs_.forcey;
		ghost = lhs_.ghost;
		id = lhs_.id;
		Jxy = lhs_.Jxy;
		Jyx = lhs_.Jyx;
		Jxx = lhs_.Jxx;
		Jyy = lhs_.Jyy;
		
		// A copy operation does not copy the trajectories. The trajectories
		// are only created for new particles.
		
		/*numTrajPoints=lhs_.numTrajPoints;
		
		for (std::vector<double>::iterator p=(*(lhs_.trajectoryx)).begin();
				p!=(*(lhs_.trajectoryx)).end(); ++p)
		{
			(*trajectoryx).push_back(*p);
		}
		
		for (std::vector<double>::iterator p=(*(lhs_.trajectoryy)).begin();
				p!=(*(lhs_.trajectoryy)).end(); ++p)
		{
			(*trajectoryy).push_back(*p);
		}*/
		
		//return *this;
				
	//};
	
	void set_Tracked();								// Sets particle tracking for this particle
	
	bool get_Tracked() const;					// Returns Tracked
	
	double get_force_dx() const { return force_dx;}
	
	double get_force_dy() const { return force_dy;}
	
	double get_force_tx() const { return force_tx;}
	
	double get_force_ty() const { return force_ty;}
	
	void set_force_d_t( double force_dx_, double force_dy_, double force_tx_, double force_ty_ )
	{
		force_dx = force_dx_;
		force_dy = force_dy_;
		force_tx = force_tx_;
		force_ty = force_ty_;
	}
	
	
};



#endif

