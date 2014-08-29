#ifndef CCOORD_HPP
#define CCOORD_HPP


class CCoord {
 public:
    double x, y;
    
  public:
	CCoord() {
		x=0;
		y=0;
				
	}
	
	~CCoord() {}
	
	CCoord operator+ (CCoord param) {
		CCoord temp;
		temp.x = x + param.x;
		temp.y = y + param.y;
		return (temp);
	};
	
	CCoord operator- (CCoord param) {
		CCoord temp;
		temp.x = x - param.x;
		temp.y = y - param.y;
		return (temp);
	};
	
	/*bool operator <(CCoord const a) const
	{
		return x < a.x && y < a.y;
		
	//	return x < a.x;
	};
	
	bool operator ==(CCoord const a) const
	{
		return x == a.x && y == a.y;
		
	//	return x < a.x;
	};
	
	bool operator !=(CCoord const a) const
	{
		return x != a.x || y != a.y;
		
	//	return x < a.x;
	};
	
	bool operator >(CCoord const a) const
	{
		return x > a.x;
		
	//	return x < a.x;
	};
	
	
	*/
	void set_coords (double a,double b)	{
			x = a;
			y = b;
			
	};
    
    double get_x ()	{
			return x;
	};
    
    double get_y ()	{
			return y;
	}
    
    
};


#endif
