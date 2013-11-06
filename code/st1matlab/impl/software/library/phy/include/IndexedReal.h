//
//  IndexedReal.h
//  Utility class for storing a real and an integer, typically it's index.
//  This is used for sorting and keeping track of permutations.
//
//  Created by Robby McKilliam on 5/06/12.
//  Copyright (c) 2012 ITR, University of South Australia. All rights reserved.
//

#ifndef _IndexedReal_h_included_
#define _IndexedReal_h_included_

#ifndef MSGTYPE
#define MSGTYPE double
//#define MSGTYPE float
#endif

class IndexedReal
{

public:
  IndexedReal( MSGTYPE value, int index );

  bool operator<(const IndexedReal& other) const;
  
  MSGTYPE v;
  int i;

private:
  int sigma(int k);

};

#endif
