#ifndef CDELTRIANGLE_HPP
#define CDELTRIANGLE_HPP

class CDelTriangle {
    double Ax,Ay,Bx,By,Cx,Cy; 
    int daughter1,daughter2, daughter3;
    bool divided;    
    bool finalDaughter;
    
  public:
	CDelTriangle() {
		Ax=0;
		Ay=0;
		Bx=0;
		By=0;
		Cx=0;
		Cy=0;
		divided=false;
		daughter1=-1;
		daughter2=-1;
		daughter3=-1;
		finalDaughter=true;
	}
	
	~CDelTriangle() {}
	
		 
    void set_vertices (double valAx, double valAy, double valBx, double valBy, double valCx, double valCy)	{
			Ax=valAx;
			Ay=valAy;
			Bx=valBx;
			By=valBy;
			Cx=valCx;
			Cy=valCy;
			
	};
	
	void set_divided() {
			divided=true;
		
	};
	
	void set_finalDaughter(bool val) {
			finalDaughter=val;
		
	};
	
	void set_daughters(int valdaughter1, int valdaughter2, int valdaughter3) {
		daughter1 = valdaughter1;
		daughter2 = valdaughter2;
		daughter3 = valdaughter3;
		finalDaughter=false;
	};
	
	int get_daughter1 ()	{
			return daughter1;
	};
	
	int get_daughter2 ()	{
			return daughter2;
	};
	
	int get_daughter3 ()	{
			return daughter3;
	};
	
	bool get_finalDaughter ()	{
			return finalDaughter;
	};
	
	bool get_divided ()	{
			return divided;
	};
	
		double get_Ax ()	{
			return Ax;
	};
    
    	double get_Ay ()	{
			return Ay;
	};
	
		double get_Bx ()	{
			return Bx;
	};
	
		double get_By ()	{
			return By;
	};
	
	double get_Cx ()	{
	
			return Cx;
	};
	double get_Cy ()	{
			return Cy;
	}
    
  
    
};
#endif
