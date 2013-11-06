% CoherentMackenthun   MATLAB wrapper for C++ implementation of the
%    Coherent Mackenthun phase estimator. 
%
%   Operate as follows:  Firstly construct an estimator with
%
%   CoherentMackenthun('name',D,P,p)
%   
%   where D, and P and vectors of positive integers representing
%   the position of the pilot and data symbols, p is a vector of
%   complex number representing the pilot symbols.  'name' can be
%   any string.  This function allocates memory for the estimator.  
%   Calling it again will reallocate memory.  
%
%   Run the estimator with
%
%   ahat = CoherentMackenthun('name', y)
%   
%   where y is a vector of recieved complex symbols and ahat is the
%   estiamted complex gain.  y must be of length |P|+|D|.
%
%   Multiple estimators with different names can coexist.
%
%   AUTHOR: Robby McKilliam (robby.mckilliam@unisa.edu.au)
%   Copyright 2010-2012 Robby McKilliam ITR-UniSA
%