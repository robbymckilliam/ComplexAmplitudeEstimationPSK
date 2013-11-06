#include "IndexedRealSDR.h"

/// \brief	Standard constructor.
///
IndexedReal::IndexedReal( CMACK_MSGTYPE value, int index ) : v(value), i(index)
{
}

/// Overload the < operator for sorting.
bool IndexedReal::operator<(const IndexedReal& other) const 
{
  return v < other.v;
}


