#ifndef CLINEIDS_HPP
#define CLINEIDS_HPP

class CLineIDs
{
public:
	int id1, id2;
	
	CLineIDs(){};
	
	CLineIDs(int id1_, int id2_)
	{
		if(id1_<id2_)
		{
			id1=id1_;
			id2=id2_;
		}
		else
		{
			id2=id1_;
			id1=id2_;
		}	
		
	};
	
	~CLineIDs(){};
	
	bool operator< (const CLineIDs& lhs_)
	{
		if(lhs_.id1 == id1)
				return lhs_.id2 >= id2;
		else if(lhs_.id1 != id1)
				return lhs_.id1 >= id1;
			
		return false;
	};
	
	friend bool operator< (const CLineIDs& lhs_, const CLineIDs& rhs_)
	{
		if(lhs_.id1 == rhs_.id1)
				return lhs_.id2 >= rhs_.id2;
		else if(lhs_.id1 != rhs_.id1)
				return lhs_.id1 >= rhs_.id1;
			
		return false;
	};
	
	bool operator== (const CLineIDs& lhs_)
	{
		return (lhs_.id1 == id1 && lhs_.id2 == id2);
				
	};
	
	friend bool operator== (const CLineIDs& lhs_, const CLineIDs& rhs_)
	{
		return (lhs_.id1 == rhs_.id1 && lhs_.id2 == rhs_.id2);
				
	};
	
	
};

#endif
