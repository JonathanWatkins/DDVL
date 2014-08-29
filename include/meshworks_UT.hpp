#ifndef _MESHWORKS_UT_H_
	#define _MESHWORKS_UT_H_

using namespace std;

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

template<class T>
double gen_normal_3(T &generator)
{
  double rn=0;
  do {
		rn =generator();
	}
	while ( boost::math::isnan((double)rn)==true || boost::math::isinf((double)rn)==true  || rn >1 || rn <-1);
  return rn; 
}

boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),
              boost::normal_distribution<>());



bool eqtest( double a, double b) {
	double EQtollerance=1e-5; //was 1e-5
	if (fabs(a/b -1) < EQtollerance){
		 return true; }
	else
	{ return false; }
}


double gaussianRand() {
	 double result;
	 do {
		double u1 = (float)rand()/(float)RAND_MAX;
		double u2 = (float)rand()/(float)RAND_MAX;
		result = sqrt((double)-2*log(u1))*cos(2*pi*u2);
		}
	 while ( boost::math::isnan((double)result)==true || boost::math::isinf(result)==true  || result >1 || result <-1);
	 if (result != result) {
	cout << "ERR: gaussian rand result: " << "(" << result << ")" << endl;
	}
	if (boost::math::isinf(result)) cout << "inf: gaussian rand result: " << "(" << result << ")" << endl;
	

	//else cout << "gaussian: "  << result << endl;
	 
	 return result;
}


inline double forceForm(double dist, double f0, double lambda ) {
	double force=0;
		
	
	if (dist==0) {
		cout << "zero" << endl;
		force=-0.0000000001*fabs(gaussianRand());
	}
	else {
		force = f0*boost::math::cyl_bessel_k(1,  dist/lambda);
	}
	
	return force;
	
}

void wrapVortices(list<CVortex>& vorticesList, double forceRange, int t) {
	list<CVortex> wrappedVorticesList;
	wrappedVorticesList=vorticesList;
  for (list<CVortex>::iterator p = vorticesList.begin();
			p!=vorticesList.end(); ++p ) {
		// wrap vortices on tube
		if (p->get_x() >=0) {
			if (p->get_y() <= forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelHeight+b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= channelHeight-forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelHeight-b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
		}
		else if (p->get_x() < 0) {  // wrap vortices on cone
			if (p->get_y() <= -fabs(p->get_x())*tan(funnelAngle)+forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelHeight+b0+2.0*fabs(p->get_x())*tan(funnelAngle) );
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= channelHeight+fabs(p->get_x())*tan(funnelAngle)-forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelHeight-b0-2.0*fabs(p->get_x())*tan(funnelAngle));
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			
		}
		
		
		
	}
	vorticesList=wrappedVorticesList;
	
}



void calculateForcesSerial(list<CVortex> vorticesList, list<CVortex>& pinsList, list<CVortex>::iterator q, int t, double tau, double eta, double kB, double dt, double forceRange,double f0, double lambda, double b0) {
	
	//list<CVortex> wrappedVorticesList;
	if (tube==geometry) wrapVortices(vorticesList, forceRange, t);
	
	
	double JyyK=0;  // kinetic term
	double JyyV=0;  // potential term
	
	double JxxK=0;  // kinetic term
	double JxxV=0;  // potential term
	
	double JxyK=0;  // kinetic term
	double JxyV=0;  // potential term
	
	double JyxK=0;  // kinetic term
	double JyxV=0;  // potential term
	
	
		
	double forcex=0;
	double forcey=0;
	
	// temp solution
	double v1[2]={0,0};
	double v2[2]={0,0}; 
	double vortexSum[2]={0,0};
	double pinningSum[2]={0,0};
	double directionVector[2]={0,0};
	double tempSum[2]={0,0};
	
	
	JyyK=0;// -1*q->get_vely()*q->get_vely();
	JxxK=0;// -1*q->get_velx()*q->get_velx();
	JxyK=0;// -1*q->get_velx()*q->get_vely();
	JyxK=0;// -1*q->get_velx()*q->get_vely();
	
	
	
	// vortex sum
	for (list<CVortex>::iterator p = vorticesList.begin();
			p != vorticesList.end(); ++p) {
		
		if (q->get_id()==p->get_id())
			continue;
				
				v1[0]=q->get_x();
				v1[1]=q->get_y();
				v2[0]=p->get_x();
				v2[1]=p->get_y();
				
				double distcal= sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]));
			  
			  if (distcal!=distcal) { // error checking
					cout << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cout << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					cerr << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cerr << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					
					}
				if (boost::math::isinf(distcal))
				cout << "t: " << t << "dist inf" << "(" << distcal << ")" << endl;
		
					
			if (distcal > forceRange) //only include close vortices
				continue;
				
				directionVector[0]=v2[0]-v1[0];
				directionVector[1]=v2[1]-v1[1];
			  double mod = sqrt((double)directionVector[0]*directionVector[0]+directionVector[1]*directionVector[1]);
				
				if (0==mod) {
					
					vortexSum[0]=vortexSum[0]-0.0000000001*(rand()%3-1);
					vortexSum[1]=vortexSum[1]-0.0000000001*(rand()%3-1);
					cout << "Two vortices at identical positions" << endl;
				}
				else {
					directionVector[0]=directionVector[0]/mod;
					directionVector[1]=directionVector[1]/mod;
					double f=forceForm(distcal,f0, lambda);
					vortexSum[0]=vortexSum[0]+f*directionVector[0];
					vortexSum[1]=vortexSum[1]+f*directionVector[1];
					
					JyyV+=0.5*(p->get_y()-q->get_y())*f*directionVector[1];
					JxxV+=0.5*(p->get_x()-q->get_x())*f*directionVector[0];
					JxyV+=0.5*(p->get_y()-q->get_y())*f*directionVector[0];
					JxyV+=0.5*(p->get_x()-q->get_x())*f*directionVector[1];
				
				
				}
				
				if (vortexSum[0]!=vortexSum[0] || vortexSum[1] != vortexSum[1])
				cout << "t: " << t << "partial sum vort nan" << "(" << vortexSum[0] << ", " << vortexSum[1] << ") mod: " << mod << endl;
				if (boost::math::isinf(vortexSum[0]) || boost::math::isinf(vortexSum[1]))
				cout << "t: " << t << "partial sum vort inf" << "(" << vortexSum[0] << ", " << vortexSum[1] << ") mod: " << mod << endl;
	
				
				
			
			
		}
			
		
		
		
	if (vortexSum[0]!=vortexSum[0] || vortexSum[1] != vortexSum[1])
	cout << "t: " << t << "vort nan" << "(" << vortexSum[0] << ", " << vortexSum[1] << ")" << endl;
	if (boost::math::isinf(vortexSum[0]) || boost::math::isinf(vortexSum[1]))
	cout << "t: " << t << "vort inf" << "(" << vortexSum[0] << ", " << vortexSum[1] << ")" << endl;
	
	
	
	// pin sum
	
	for (list<CVortex>::iterator p = pinsList.begin();
			p != pinsList.end(); ++p) {
	      v1[0]=q->get_x();
				v1[1]=q->get_y();
				v2[0]=p->get_x();
				v2[1]=p->get_y();
			
			double distcal= sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]));
			
			if (distcal!=distcal) { // error checking
					cout << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cout << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					cerr << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cerr << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					
			}
			
			if (distcal > forceRange) //only include close vortices
				continue;
					
					
				directionVector[0]=v2[0]-v1[0];
				directionVector[1]=v2[1]-v1[1];
				double mod = sqrt((double)directionVector[0]*directionVector[0]+directionVector[1]*directionVector[1]);
				if (0!=mod) {
					directionVector[0]=directionVector[0]/mod;
					directionVector[1]=directionVector[1]/mod;
				}
				
				double f=forceForm(distcal,f0, lambda);
					
				pinningSum[0]=pinningSum[0]+f*directionVector[0];
				pinningSum[1]=pinningSum[1]+f*directionVector[1];
				
				JyyV+=0.5*(p->get_y()-q->get_y())*f*directionVector[1];
			  JxxV+=0.5*(p->get_x()-q->get_x())*f*directionVector[0];
				JxyV+=0.5*(p->get_y()-q->get_y())*f*directionVector[0];
				JxyV+=0.5*(p->get_x()-q->get_x())*f*directionVector[1];
				
			
		}
	if (pinningSum[0]!=pinningSum[0] || pinningSum[1] != pinningSum[1])
	cout << "t: " << t << "pin nan" << "(" << pinningSum[0] << ", " << pinningSum[1] << ")" << endl;
	if (boost::math::isinf(pinningSum[0]) || boost::math::isinf(pinningSum[1]))
	cout << "t: " << t << "pin inf" << "(" << pinningSum[0] << ", " << pinningSum[1] << ")" << endl;
	
	
	
	
		tempSum[0]=sqrt((double)tau)*gen_normal_3(generator);
		tempSum[1]=sqrt((double)tau)*gen_normal_3(generator);
		
	
	
	if (tempSum[0]!=tempSum[0] || tempSum[1] != tempSum[1]) {
	  cout << "t: " << t << "temp nan" << "(" << tempSum[0] << ", " << tempSum[1] << ")" << endl;
		tempSum[0]=0;
	  tempSum[1]=0;
	}
	if (boost::math::isinf(tempSum[0]) || boost::math::isinf(tempSum[1]))
	cout << "t: " << t << "temperature inf" << "(" << tempSum[0] << ", " << tempSum[1] << ")" << endl;
	
	
	
	forcex= vortexSum[0]+pinningSum[0];
	forcey= vortexSum[1]+pinningSum[1];
	
	 
	if (forcex!=forcex || forcey != forcey || tempSum[0] != tempSum[0]
	 || tempSum[1]!=tempSum[1])
	cout << "t: " << t << "force nan" << "(" << forcex << ", " << forcey << ")" << endl;
	if (boost::math::isinf(forcex) || boost::math::isinf(forcey))
	cout << "t: " << t << "force inf" << "(" << forcex << ", " << forcey << ")" << endl;
	
		
	
	double velx=(forcex+tempSum[0])/eta;
	double vely=(forcey+tempSum[1])/eta;
	if (q->get_x() >sourceWidth+channelWidth) { // rescaled viscosity for sink vortices
		velx=velx*2;
		vely=vely*2;
	}
	if (velx!=velx) cout << "t: " << t << " velx nan" << endl;
	if (vely!=vely) cout << "t: " << t << " vely nan" << endl;
	if (boost::math::isinf(velx) || boost::math::isinf(vely)) cout << "t: " << t << " vel nan " << velx << ", " << vely << endl; 
	q->set_vel(velx,vely);
	//cout << "dt=" << dt << endl;
	//cout << "v passed to object (" << velx << ", " << vely << ")" << endl; 
	//cout << "MAX allowed dr/dt: " << 0.5*b0/dt<< endl;
	//cout << "v returned from object (" << q->get_velx() << ", " << q->get_vely() << ")" << endl; 
	
	// Velocities are allowed to be a maximum 0.5*b0 /dt to avoid escaping vortices
	
	if(q->get_velx()>0.5*b0/dt) {
		//cout << "  vx velocity rectified " << q->get_velx(); 
		q->set_velx(0.5*b0/dt);
		//cout << " changed to " << setw(15) << q->get_velx() << endl;
	}
	
	if(q->get_velx()<-0.5*b0/dt) {
		//cout << "  vx velocity rectified " << q->get_velx(); 
		q->set_velx(-0.5*b0/dt);
		//cout << " changed to " << setw(15) << q->get_velx() << endl;
	}
	
	
	if(q->get_vely()>0.5*b0/dt) {
		//cout << "  vy velocity rectified " << q->get_vely(); 
	  q->set_vely(0.5*b0/dt);
	  //cout << " changed to "<< setw(15) << q->get_vely() << endl;
	}
	
	if(q->get_vely()<-0.5*b0/dt) {
		//cout << "  vy velocity rectified " << q->get_vely(); 
	  q->set_vely(-0.5*b0/dt);
	  //cout << " changed to "<< setw(15) << q->get_vely() << endl;
	}
	
	
	
	//cout << "pos: " << q->get_x() << ", " << q->get_y() << endl;
	q->set_pos(q->get_x()+q->get_velx()*dt,q->get_y()+q->get_vely()*dt);
	
	//if (fabs(q->get_velx()*dt)>0.02*a0) cout << "t: " << t << " drx=" << 100*q->get_velx()*dt/a0 << "% of a0" << endl;  
	
	//if (fabs(q->get_vely()*dt)>0.02*a0) cout << "t: " << t << " dry=" << 100*q->get_vely()*dt/a0 << "% of a0" << endl;  
	
	if (q->get_velx()*dt>0.5*b0 || q->get_velx()*dt<-0.5*b0 ) {
		cout << "vibration bigger than 0.5b0 in a single timestep!" << endl;
	}
  
  if (q->get_vely()*dt>0.5*b0 || q->get_vely()*dt<-0.5*b0 ) {
		cout << "vibration bigger than 0.5b0 in a single timestep!" << endl;
	}
  
  // Set partilces local stress
  
  q->set_Jyy(JyyK+JyyV);
  q->set_Jxx(JxxK+JxxV);
  q->set_Jxy(JxyK+JxyV);
  q->set_Jyx(JyxK+JyxV);
  
  
   
  
	
}







#endif
