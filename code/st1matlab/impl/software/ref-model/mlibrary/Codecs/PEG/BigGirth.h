#include <stdlib.h>
#include <iostream> // C++ I/O library header

#ifndef BIGGIRTH
#define BIGGIRTH
class Random{
private:
	
  unsigned long int seed;  //previously LONG INT
  unsigned long int seed_u;
  
public:
	
  Random(void){srand((unsigned int)time(NULL));}
  ~Random(void){;}
  long int uniform(long int a, long int b);	
}; 

class Node
{
public:
	int	d;
	int dmax;
	int	*connection;
	
	Node(void);
	~Node(void);
	void Init(int deg, int max);
	void Init(int deg);
	
	bool active();
};


class BigGirth
{
	public:
		long int						M, N;
		long int						K;
		long int						EXPAND_DEPTH;
		char					*filename;
		long int						*localGirth;
		Node					*varnodes, *chknodes;
		Random				*myrandom;

		BigGirth(long int m, long int n, int *symbolDegSequence, int *checkDegSequence, char *filename, long int sglConcent, long int tgtGirth);
		BigGirth(void);

		void	writeToFile_Hcompressed(void);

		~BigGirth(void);
	
	private:
		long int		selectParityConnect(int kthSymbol, int mthConnection, long int &cycle);
		void				updateConnection(int kthSymbol);
};
#endif





