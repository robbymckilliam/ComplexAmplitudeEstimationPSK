#include "IndexedReal.h"

/// \brief	Standard constructor.
///
IndexedReal::IndexedReal( MSGTYPE value, int index ) : v(value), i(index)
{
}

/// Overload the < operator for sorting.
bool IndexedReal::operator<(const IndexedReal& other) const 
{
  return v < other.v;
}


