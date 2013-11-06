//
//  IndexedReal.h
//  Utility class for storing a real and an integer, typically it's index.
//  This is used for sorting and keeping track of permutations.
//
//  Created by Robby McKilliam on 5/06/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _IndexedRealSDR_h_included_
#define _IndexedRealSDR_h_included_

#ifndef CMACK_MSGTYPE
#define CMACK_MSGTYPE double
//#define MSGTYPE float
#endif

class IndexedReal
{

public:
  IndexedReal( CMACK_MSGTYPE value, int index );

  bool operator<(const IndexedReal& other) const;
  
  int i;
  CMACK_MSGTYPE v;

private:
  int sigma(int k);

};

#endif
