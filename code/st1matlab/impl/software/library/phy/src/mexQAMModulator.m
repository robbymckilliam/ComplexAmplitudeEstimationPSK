% mexQAMModulator   MATLAB wrapper for QAMModulator class
%    This mex-file enables the use of the QAMModulator class from MATLAB
%
%    Syntax:
%    The first argument is the ID of the QAM Modulator (starting from 0). This
%    allows to use many modulators in parallel. The maximum number of
%    parallel modulators is defined as MAXMOD (default=100).
%
%    initialise:
%       mexQAMModulator(ID, 'init', length, ModulationMap)
%
%         length        ... number of transmit symbols
%         ModulationMap ... vector with complex constellation points
%
%    set channel probabilities:
%       mexQAMModulator(ID, 'set-ch', y, sigma2)
%
%         y             ... complex received symbol of length "length"
%         sigma2        ... variance of AWGN
%
%    set a-priori probabilities:
%       mexQAMModulator(ID, 'set-apriori', La)
%
%       La is a matrix with "bitsPerSymbol" rows and "length" columns holding LLRs
%
%    modulate:
%       x = mexQAMModulator(ID, 'modulate', bits)
%
%       Maps bits to complex transmit symbols.
%
%    demodulate:
%       Le = mexQAMModulator(ID, 'demodulate')
%
%       Returns a matrix with "bitsPerSymbol" rows and "length" columns holding LLRs
%       based on the previously set channel and a-priori messages.
%
%    remodulate:
%       [xhat, resvar] = mexQAMModulator(ID, 'remodulate')
%
%       Returns a vector of length "length" with complex remodulated
%       symbols and the residual variance.
%
%    modulator configuration:
%       mexQAMModulator(ID, 'configure', parameter, value)
%
%       Possible parameters:
%
%   (C) 2012 Gottfried Lechner, ITR/UniSA
%