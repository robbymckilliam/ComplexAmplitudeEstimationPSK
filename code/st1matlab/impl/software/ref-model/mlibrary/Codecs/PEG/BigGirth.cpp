#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "BigGirth.h"

#define MAX(a,b)	((a)>(b) ? (a) : (b))
#define MIN(a,b)	((a)<(b) ? (a) : (b))

using namespace std;

long int Random::uniform(long int  a, long int b) // [a, b-1]
{
  return (long int)((double)rand()/(RAND_MAX+1.0)*(b-a)+a);
}

Node::Node(void)
{
	connection = NULL;
}

Node::~Node(void)
{
	delete connection;
}

void Node::Init(int deg, int max)
{
	d						= deg;
	dmax				= max;
	connection	= new int[dmax];
}

void Node::Init(int deg)
{
	Init(deg, deg);
}

bool Node::active()
{
	return d<dmax;
}


// ============================================================================

BigGirth::BigGirth(void)
{
}

BigGirth::BigGirth(long int rows, long int cols, int *symbolDegSequence, int *checkDegSequence, char *filename, long int sglConcent, long int tgtGirth)
{
	long int n, m, d, index, localDepth=2*cols;
	int iter;

	EXPAND_DEPTH			= MAX((tgtGirth-4)/2, 0);
	myrandom					= new Random();
	M									= rows;
	N									= cols;
	(*this).filename	= filename;
	localGirth				= new long int[N];
	varnodes					= new Node[N];
	chknodes					= new Node[M];

	for(n=0;n<N;n++)
		varnodes[n].Init(symbolDegSequence[n]);
	
	for(m=0;m<M;m++)
	{
		if(sglConcent==0)
			chknodes[m].Init(0, checkDegSequence[m]);
		else
			chknodes[m].Init(0, 2*checkDegSequence[m]);
	}
      
	for(n=0;n<N;n++)
	{
		index	= -1;
		for(m=0;m<M;m++)
		{
			if(((chknodes[m].d<d) || (index<0)) && chknodes[m].active())
			{
				d			= chknodes[m].d;
				index	= m;
			}
		}
		varnodes[n].connection[0] = index;

		iter = 0; 
		ITER:
		localGirth[n] = 2*N;
		for(d=1;d<varnodes[n].d;d++)
		{
			varnodes[n].connection[d] = selectParityConnect(n, d, localDepth);
			localGirth[n] = MIN(localGirth[n], localDepth);
			
			if(n>0 && localGirth[n]<localGirth[n-1] && iter<20)
			{
				iter++;
				goto ITER;
			}
			if(localGirth[n]==0 && iter<30)
			{
				iter++;
				goto ITER;
			}
		}
		if((n%100)==0)
		{
			cout<<"n="<<n<<" \tdeg="<<varnodes[n].d<<" \tLocalGirth=";
			if(localGirth[n]==2*N)
				cout<<"Inf";
			else
				cout<<2*localGirth[n]+4;
			cout<<endl;
		}

		updateConnection(n);
	}

	cout<<"Showing the row weight distribution..."<<endl;
	for(m=0;m<M;m++)
		cout<<chknodes[m].d<<" ";
	cout<<endl;

  ofstream	cycleFile;
	char			lhgName[100];
	
	sprintf(lhgName, "%s.lhg", filename);
	cycleFile.open(lhgName,ios::out);
	localDepth = 2*N;
  for(n=0;n<N;n++)
	{
    if(localGirth[n]<localDepth)
			localDepth = localGirth[n];
		if(localDepth==2*N)
			cycleFile<<"Inf ";
		else
			cycleFile<<2*localDepth+4<<" ";
  }
  cycleFile<<endl;
  cycleFile.close();	
	
	cout<<"*************************************************************"<<endl;
	cout<<"       The global girth of the PEG Tanner graph :="<< 2*localDepth+4<<endl;
	cout<<"*************************************************************"<<endl;
}

BigGirth::~BigGirth(void)
{
	delete localGirth;
	delete myrandom;
}

long int BigGirth::selectParityConnect(int kthSymbol, int mthConnection, long int &cycle)
{
	long int i, j, k, index, mincycles, numCur, cpNumCur;
	long int *tmp;
	long int *current;

	mincycles	= 0;
	tmp				= new long int[M];
  
	numCur		= mthConnection;
	current		= new long int[mthConnection];
	for(i=0;i<mthConnection;i++)
		current[i] = varnodes[kthSymbol].connection[i];

	LOOP:
	mincycles++;
	for(i=0;i<M;i++)
		tmp[i] = 0; 
	for(i=0;i<mthConnection;i++)
		tmp[varnodes[kthSymbol].connection[i]] = 1;
	for(i=0;i<numCur;i++)
		for (j=0; j<chknodes[current[i]].d; j++)
			for(k=0;k<varnodes[chknodes[current[i]].connection[j]].d;k++)
				tmp[varnodes[chknodes[current[i]].connection[j]].connection[k]] = 1;

	index			= 0;
	cpNumCur	= 0;
	for(i=0;i<M;i++)
	{
		if(tmp[i]==1)
			cpNumCur++;
		if(tmp[i]==1 || chknodes[i].d>=chknodes[i].dmax) 
			index++;   
	}

	if(cpNumCur==numCur)
	{
		//can not expand any more
		//additional handlement to select one having least connections
		j=10000000; //dummy number
		for(i=0;i<M;i++)
			if(tmp[i]==0 && chknodes[i].d<j && chknodes[i].active())
				j = chknodes[i].d;

		for(i=0;i<M;i++)
			if(tmp[i]==0)
				if(chknodes[i].d!=j || !chknodes[i].active())
					tmp[i] = 1;

		index = 0;
		for(i=0;i<M;i++)
			index += tmp[i];
    
		j			= myrandom->uniform(0, M-index)+1;
		index	= 0;
		for(i=0;i<M;i++)
		{
			if(tmp[i]==0)
				index++;
			if(index==j)
				break;
		}
		delete tmp;
		tmp			= NULL;
		delete current;
		current	= NULL;
		return i;
	}
	else
	{
		if(index==M || mincycles>EXPAND_DEPTH)
		{
			//covering all parity nodes or meet the upper bound on cycles
			cycle = mincycles-1;
			for(i=0;i<M;i++)
				tmp[i]					= 0;
			for(i=0;i<numCur;i++)
				tmp[current[i]]	= 1;

			index=0;
			for(i=0;i<M;i++)
				if(tmp[i]==1)
					index++;
			if(index!=numCur)
			{
				cout<<"Error in the case of (index==M)"<<endl;
				exit(-1);
			}
    
			//additional handlement to select one having least connections
			j=10000000; 
			for(i=0;i<M;i++)
				if(tmp[i]==0 && chknodes[i].d<j && chknodes[i].d<chknodes[i].dmax)
					j = chknodes[i].d;

			for(i=0;i<M;i++)
				if(tmp[i]==0)
					if(chknodes[i].d!=j || chknodes[i].d>=chknodes[i].dmax)
						tmp[i] = 1;
      
			index = 0;
			for(i=0;i<M;i++)
				if(tmp[i]==1)
					index++;
   
			j			= myrandom->uniform(0, M-index)+1;
			index = 0;
			for(i=0;i<M;i++)
			{
				if(tmp[i]==0)
					index++;
				if(index==j)
					break;
			}
			delete tmp;
			tmp			= NULL;
			delete current;
			current	= NULL;
			return i;
		}
		else
		{
			if(cpNumCur>numCur && index!=M)
			{
				delete current;
				current	= NULL;
				numCur	= cpNumCur;
				current	= new long int[numCur];
				index		= 0;
				for(i=0;i<M;i++)
					if(tmp[i]==1)
					{
						current[index]=i;
						index++;
					}
				goto LOOP;
			}
			else 
			{
				cout<<"Should not come to this point..."<<endl;
				cout<<"Error in BigGirth::selectParityConnect()"<<endl;
				return(-1);
			}
		}
	}
}


void BigGirth::updateConnection(int kthSymbol)
{
	int i, m;

	for(i=0;i<varnodes[kthSymbol].d;i++)
	{
		m		= varnodes[kthSymbol].connection[i];
		chknodes[m].connection[chknodes[m].d] = kthSymbol;
		chknodes[m].d++;
	}
}

void BigGirth::writeToFile_Hcompressed(void){
  int i, j, max_col;
  int *(*parityCheck_compressed);

  //cout<<"---------------code format--------------------------"<<endl;
  //cout<<"-            Block length N                        -"<<endl;
  //cout<<"-            Num of Check Nodex M                  -"<<endl;
  //cout<<"-            Num of column in the compressed H     -"<<endl;
  //cout<<"-            H matrix (compressed)                 -"<<endl;
  //cout<<"----------------------------------------------------"<<endl;

  //finding the num of columns, l, of the compressed parity-check matrix	
  max_col = 0;
  for(i=0;i<M;i++)
    if(chknodes[i].d>max_col) 
      max_col = chknodes[i].d;

  parityCheck_compressed	= new int* [M];
  for(i=0;i<M;i++)
    parityCheck_compressed[i] = new int[max_col];
  for(i=0;i<M;i++)
	{
    for(j=0;j<max_col;j++)
			parityCheck_compressed[i][j] = 0;
    for(j=0;j<chknodes[i].d;j++)
	    parityCheck_compressed[i][j] = chknodes[i].connection[j]+1; 
  }

  ofstream codefile;  
  codefile.open(filename,ios::out);
  codefile<<N<<endl;
  codefile<<M<<endl;
  codefile<<max_col<<endl;
  for(i=0;i<M;i++)
	{
    for(j=0;j<max_col;j++)
      codefile<<parityCheck_compressed[i][j]<<" ";
    codefile<<endl;
  }
  codefile.close();

  for(i=0;i<M;i++)
    delete parityCheck_compressed[i];
  delete parityCheck_compressed;
}


