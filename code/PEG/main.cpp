/************************************************************************/
/*                                                                      */
/*        Free software: Progressive edge-growth (PEG) algorithm        */
/*        Created by Xiaoyu Hu                                          */
/*                   Evangelos Eletheriou                               */
/*                   Dieter Arnold                                      */
/*        IBM Research, Zurich Research Lab., Switzerland               */
/*                                                                      */
/*        The C++ sources files have been compiled using xlC compiler   */
/*        at IBM RS/6000 running AIX. For other compilers and platforms,*/
/*        minor changes might be needed.                                */
/*                                                                      */
/*        Bug reporting to: xhu@zurich.ibm.com                          */
/************************************************************************/

/************************************************************************/
/*                                                                      */
/*        Modifications made by Gottfried Lechner                       */
/*                                                                      */
/*        - fixed a memory leak that lead to a crash for long codes     */
/*        - added the option to supply the check node distribution      */
/*        - changes such that codes with M>N can be constructed         */
/*        - restructuring of the code                                   */
/*        - changed parameter checking to ignore case                   */
/*                                                                      */
/*        Bug reporting to: gottfried.lechner@unisa.edu.au              */
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "BigGirth.h"

using namespace std;

#define EPS  1e-6

int main(int argc, char * argv[]){
  long int	i, j, m, mc, N, M, n, ev, ec;
  int				sglConcent		= 1;								// default to non-strictly concentrated parity-check distribution
  long int	targetGirth;
  char			codeName[100], degFileName[100], degFileNameChk[100];
  int				*degSeq, *deg, *degSeqChk;
  double		*degFrac;
  BigGirth	*bigGirth;
	
  int numArgs=(argc-1)/2;
  if (argc<9) {
  USE:
    cout<<"*******************************************************************************************"<<endl;
    cout<<" Usage Reminder: MainPEG -numM M -numN N -codeName CodeName -degFileName DegFileName       "<<endl;
    cout<<"         option:         -sglConcent SglConcent -degFileNameChk ChkDegFileName             "<<endl; 
    cout<<"                         sglConcent==0 ----- strictly concentrated parity-check            "<<endl;
    cout<<"                                       degree distribution (including regular graphs)      "<<endl;
    cout<<"                         sglConcent==1 ----- Best-effort concentrated (DEFAULT)            "<<endl;
    cout<<"         option:         -tgtGirth TgtGirth                                                "<<endl; 
    cout<<"                  TgtGirth==4, 6 ...; if very large, then greedy PEG (DEFAULT)             "<<endl;
    cout<<"                  IF sglConcent==0, TgtGirth is recommended to be set relatively small     "<<endl;
    cout<<"                                                                                           "<<endl;
    cout<<" Remarks: File CodeName stores the generated PEG Tanner graph. The first line contains     "<<endl;
    cout<<"          the block length, N. The second line defines the number of parity-checks, M.     "<<endl;
    cout<<"          The third line defines the number of columns of the compressed parity-check      "<<endl;
    cout<<"          matrix. The following M lines are then the compressed parity-check matrix.       "<<endl;
    cout<<"          Each of the M rows contains the indices (1 ... N) of 1's in the compressed       "<<endl;
    cout<<"          row of parity-check matrix. If not all column entries are used, the column       "<<endl;
    cout<<"          is filled up with 0's.                                                           "<<endl;
    cout<<"                                                                                           "<<endl;
    cout<<"          File DegFileName is the input file to specify the degree distribution (node      "<<endl;
    cout<<"          perspective). The first line contains the number of various degrees. The second  "<<endl;
    cout<<"          defines the row vector of degree sequence in the increasing order. The vector    "<<endl;
    cout<<"          of fractions of the corresponding degree is defined in the last line.            "<<endl;
    cout<<"                                                                                           "<<endl;
    cout<<"          A log file called 'leftHandGirth.dat' will also be generated and stored in the   "<<endl;
    cout<<"          current directory, which gives the girth of the left-hand subgraph of j, where   "<<endl;
    cout<<"          1<=j<=N. The left-hand subgraph of j is defined as all the edges emanating from  "<<endl;
    cout<<"          bit nodes {1 ... j} and their associated nodes.                                  "<<endl; 
    cout<<"                                                                                           "<<endl;
    cout<<"          The last point is, when strictly concentrated parity-check degree distribution   "<<endl;
    cout<<"          is invoked, i.e. sglConcent==0, the girth might be weaken to some extent as      "<<endl;
    cout<<"          compared to the generic PEG algorithm.                                           "<<endl;
    cout<<"*******************************************************************************************"<<endl;
    exit(-1);
  }else {
    for(i=0;i<numArgs;i++){
      if (strcasecmp(argv[2*i+1], "-numM")==0) {
				M=atoi(argv[2*i+2]);
      } else if(strcasecmp(argv[2*i+1], "-numN")==0) {
				N=atoi(argv[2*i+2]);
				targetGirth = 2*N;
      } else if(strcasecmp(argv[2*i+1], "-codeName")==0) {
				strcpy(codeName, argv[2*i+2]); 
      } else if(strcasecmp(argv[2*i+1], "-degFileName")==0) {
				strcpy(degFileName, argv[2*i+2]); 
      } else if(strcasecmp(argv[2*i+1], "-degFileNameChk")==0) {
				strcpy(degFileNameChk, argv[2*i+2]); 
      } else if(strcasecmp(argv[2*i+1], "-sglConcent")==0) {
				sglConcent=atoi(argv[2*i+2]);
      } else if(strcasecmp(argv[2*i+1], "-tgtGirth")==0) {
				targetGirth=atoi(argv[2*i+2]);
      } else{
				goto USE;
      }
    }
    if(M>N) {
			cout<<"Warning: M is larger then N"<<endl;
    }
  }
	
  degSeq		= new int[N];
	degSeqChk = new int[M];
	
	// read variable node distribution
  ifstream infn(degFileName);
  if (!infn) {cout << "\nCannot open file " << degFileName << endl; exit(-1);} 
  infn >>m;
  deg			= new int[m];
  degFrac = new double[m];
  for(i=0;i<m;i++)
		infn>>deg[i];
  for(i=0;i<m;i++)
		infn>>degFrac[i];
  infn.close();  
	
	// check if it sums up to one
  double dtmp	= 0.0;
  for(i=0;i<m;i++)
		dtmp += degFrac[i];
  cout.setf(ios::fixed, ios::floatfield);
  if(fabs(dtmp-1.0)>EPS)
	{
    cout.setf(ios::fixed, ios::floatfield);
    cout <<"\n Invalid degree distribution (node perspective): sum != 1.0 but "<<setprecision(10)<<dtmp<<endl; exit(-1); 
  }
	
	// generate degree sequence for variable nodes
  for(i=1;i<m;i++)
		degFrac[i] += degFrac[i-1];
  for(i=0;i<N;i++)
	{
    dtmp	= (double)i/N;
    for(j=m-1;j>=0;j--)
		{
      if(dtmp>degFrac[j])
				break;
    }
    if(dtmp<degFrac[0])
			degSeq[i] = deg[0];
    else
			degSeq[i]	= deg[j+1];
  }
	delete deg;
  delete degFrac;

  // count number of edges induced by the variable node degree distribution
  ev=0;
  for(n=0;n<N;n++)
    ev+=degSeq[n];
  
  cout << "\nGraph has " << ev << " edges based on variable node degrees" << endl;
  
	// try to read check node distribution
  ifstream infnc(degFileNameChk);
  if (!infnc)
	{
		cout << "\nNo check node distribution provided" << endl;
    cout << "Building check node concentrated graph" << endl;
    n = ev/M;
		// generate degree sequence for check nodes
		for(m=0;m<ev-n*M;m++)
			degSeqChk[m] = n+1;
		for(m=ev-n*M;m<M;m++)
			degSeqChk[m] = n;		
	}
	else
	{
		// read check node distribution
		infnc >>mc;
		deg			= new int[mc];
		degFrac = new double[mc];
		for(i=0;i<mc;i++)
			infnc>>deg[i];
		for(i=0;i<mc;i++)
			infnc>>degFrac[i];
		infnc.close();  
		
		for(i=1;i<mc;i++)
			degFrac[i] += degFrac[i-1];
		for(i=0;i<M;i++)
		{
			dtmp	= (double)i/M;
			for(j=mc-1;j>=0;j--)
			{
				if(dtmp>degFrac[j])
					break;
			}
			if(dtmp<degFrac[0])
				degSeqChk[i] = deg[0];
			else
				degSeqChk[i] = deg[j+1];
		}
		delete deg;
		delete degFrac;
    
    // count number of edges induced by the variable node degree distribution
    ec=0;
    for(m=0;m<M;m++)
      ec+=degSeqChk[m];
    
    cout << "\nGraph has " << ec << " edges based on check node degrees" << endl;   
    
    if(ec<ev)
        cout << "\nModifying the check node degrees to match the number of edges. Check the result!" << endl;   
    
    m = M;
    while(ec<ev)
    {
      degSeqChk[--m]++;
      ec++;
    }
	}
		
	// construct code
  bigGirth	= new BigGirth(M, N, degSeq, degSeqChk, codeName, sglConcent, targetGirth);
	
	// write code
  (*bigGirth).writeToFile_Hcompressed();
	
	// cleanup
  delete degSeq;
  delete bigGirth;
}
