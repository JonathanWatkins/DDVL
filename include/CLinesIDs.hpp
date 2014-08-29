#ifndef LINESIDS_HPP
#define LINESIDS_HPP

class lineIDs
{
public:
	int id1, id2;
	
	lineIDs(){};
	
	lineIDs(int id1_, int id2_)
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
	
	~lineIDs(){};
	
	bool operator< (const lineIDs& lhs_)
	{
		if(lhs_.id1 == id1)
				return lhs_.id2 >= id2;
		else if(lhs_.id1 != id1)
				return lhs_.id1 >= id1;
	};
	
	bool operator== (const lineIDs& lhs_)
	{
		return (lhs_.id1 == id1 && lhs_.id2 == id2);
				
	};
	
	
	
};

#ENDIF
