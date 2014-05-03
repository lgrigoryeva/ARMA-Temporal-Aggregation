function A = mrdivide_my(C,B,varargin)
%MRDIVIDE Lag operator polynomial right division
%
% Syntax:
%
%   A = C/B
%	A = mrdivide(C,B,param1,val1,param2,val2,...)
%
% Description:
%
%   Given two lag operator polynomials, C(L) and B(L), perform a right
%   division so that C(L) = A(L)*B(L), or A(L) = C(L)/B(L). Right division
%   requires invertibility of the coefficient matrix associated with lag 0
%   of the denominator polynomial B(L).
%
% Input Arguments:
%
%   C - Numerator (dividend) lag operator polynomial object, as produced by
%     LagOp, in the quotient C(L)/B(L).
%
%	B - Denominator (divisor) lag operator polynomial object, as produced
%	  by LagOp, in the quotient C(L)/B(L).
%
%   If at least one of C or B is a lag operator polynomial object, the
%   other can be a cell array of matrices (initial lag operator
%   coefficients), or a single matrix (zero-degree lag operator). 
%
% Optional Input Parameter Name/Value Pairs:
%
%   'AbsTol'  Positive scalar absolute tolerance used as part of the
%             termination criterion of the calculation of the quotient
%             coefficients and, subsequently, to determine which
%             coefficients to include in the quotient. Specifying an
%             absolute tolerance allows for customization of the
%             termination criterion (see notes below). Once the algorithm
%             has terminated, 'AbsTol' is used to exclude polynomial lags
%             with near-zero coefficients. A coefficient matrix for a given
%             lag is excluded if the magnitudes of all elements of the
%             matrix are less than or equal to the absolute tolerance. The
%             default is 1e-12.
%
%   'RelTol'  Positive scalar relative tolerance used as part of the
%             termination criterion of the calculation of the quotient
%             coefficients. At each lag, a coefficient matrix is
%             calculated and its 2-norm compared to the largest coefficient
%             2-norm. If the ratio of the current norm to the largest norm
%             is less than or equal to 'RelTol', then the relative
%             termination criterion is satisfied. The default is 0.01.
%
%   'Window'  Positive integer indicating the size of the window used to
%             check termination tolerances. Window represents the number of
%             consecutive lags for which coefficients must satisfy a
%             tolerance-based termination criterion in order to terminate
%             the calculation of the quotient coefficients. If coefficients
%             remain below tolerance for the length of the specified
%             tolerance window, they are assumed to have died out
%             sufficiently to terminate the algorithm (see notes below).
%             The default is 20.
%
%   'Degree'  Nonnegative integer indicating the maximum degree of the
%             quotient polynomial. For stable denominators, the default is
%             the power to which the magnitude of the largest eigenvalue of
%             the denominator must be raised to equal the relative
%             termination tolerance 'RelTol'; for unstable denominators,
%             the default is the power to which the magnitude of the
%             largest eigenvalue must be raised to equal the largest
%             positive floating point number (see REALMAX). The default is
%             1000, regardless of the stability of the denominator.
%
% Output Argument:
%
%   A - Quotient lag operator polynomial object, with A(L) = C(L)/B(L).
%
% Notes:
%
%   o The right division operator (/) invokes MRDIVIDE, but the optional
%     inputs are available only by calling MRDIVIDE directly. 
%
%   o To right-invert a stable B(L), set C(L) = eye(B.Dimension).
%
%   o Lag operator polynomial division generally results in infinite-degree
%     polynomials. MRDIVIDE imposes a termination criterion to truncate the
%     degree of the quotient polynomial.
%
%     If 'Degree' is unspecified, the maximum degree of the quotient is
%     determined by the stability of the denominator. Stable denominator
%     polynomials usually result in quotients whose coefficients exhibit
%     geometric decay in absolute value. (When coefficients change sign, it
%     is the coefficient envelope which decays geometrically.) Unstable
%     denominators usually result in quotients whose coefficients exhibit
%     geometric growth in absolute value. In either case, maximum degree
%     will not exceed the value of 'Degree'.
%
%     To control truncation error by terminating the coefficient sequence
%     too early, the termination criterion involves three steps:
%
%      (1) At each lag in the quotient polynomial, a coefficient matrix is 
%          calculated and tested against both a relative and an absolute
%          tolerance (see 'RelTol' and 'AbsTol' inputs above).
%
%      (2) If the current coefficient matrix is below either tolerance,
%          then a tolerance window is opened to ensure that all subsequent
%          coefficients remain below tolerance for a number of lags
%          determined by 'Window'.
%
%      (3) If any subsequent coefficient matrix within the window is above
%          both tolerances, then the tolerance window is closed and
%          additional coefficients are calculated, repeating steps (1) and
%          (2) until a subsequent coefficient matrix is again below either
%          tolerance, and a new window is opened.
%
%     Steps (1)-(3) are repeated until a coefficient is below tolerance and
%     subsequent coefficients remains below tolerance for 'Window' lags, or
%     until the maximum 'Degree' is encountered, or until a coefficient
%     becomes numerically unstable (NaN or +/-Inf).
%
% References:
%
%   [1] Box, G.E.P., Jenkins, G.M., and Reinsel, G.C. Time Series Analysis: 
%       Forecasting and Control. 3rd ed. Upper Saddle River, NJ: Prentice-
%       Hall, 1994.
%
%   [2] Hayashi, F. Econometrics. Princeton, NJ: Princeton University
%       Press, 2000.
%
%   [3] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also LagOp/LagOp, LagOp/mtimes, LagOp/mldivide. 

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $   $Date: 2009/11/05 16:58:34 $

% Check input arguments:

if nargin > 10
    
   error('econ:LagOp:mrdivide:TooManyInputs',...
         'Too many inputs.')
     
end

parser = inputParser;
parser.addRequired('C');
parser.addRequired('B');
parser.addParamValue('AbsTol',LagOp.ZeroTolerance,@LagOp.checkTolerance);
parser.addParamValue('RelTol',0.01,@LagOp.checkTolerance);
parser.addParamValue('Window',20,@LagOp.checkWindow);
parser.addParamValue('Degree',[],@LagOp.checkDegree);
parser.parse(C,B,varargin{:});

C          = parser.Results.C;
B          = parser.Results.B;
absTol     = parser.Results.AbsTol;
relTol     = parser.Results.RelTol;
windowSize = parser.Results.Window;
maxDegree  = parser.Results.Degree;

if isempty(C) || isempty(B)
    
   error('econ:LagOp:mrdivide:MissingInputs',...
         'Lag operator polynomials must not be empty.')
     
end

% Allow either input to be a matrix, a cell array, or some other object, 
% and attempt a LagOp conversion. The tolerance is not applied here.

if ~isa(C,'LagOp')
   C = LagOp(C,'Tolerance',0); % Convert C(L) to a LagOp object
elseif ~isa(B,'LagOp')
   B = LagOp(B,'Tolerance',0); % Convert B(L) to a LagOp object
end

if C.Dimension ~= B.Dimension
    
   error('econ:LagOp:mrdivide:InconsistentDimension',...
         'Polynomials must have the same dimension.')
     
end

% Test the 0-lag coefficient matrix of B(L) for singularity/invertibility.
% If this matrix is noninvertible (or nearly so), then issue a single
% warning and disable the core MATLAB warnings that would otherwise be
% issued repeatedly inside the FOR loop when the "/" operator is evaluated.

wState = warning; % Save warning state

if cond(B.Coefficients{0}) > LagOp.SingularTolerance
   wString = cat(2, 'Zero-lag coefficient of denominator polynomial is \n', ...
                    '         singular or badly scaled. Results may be inaccurate.');
                
   warning('off','MATLAB:nearlySingularMatrix')
   warning('off','MATLAB:singularMatrix')
   warning('econ:LagOp:mrdivide:NearSingularMatrix',wString)
   
end

% Determine the number of quotient coefficients to calculate, and, if
% necessary, the default maximum degree of the quotient polynomial:

if isempty(maxDegree)
    
   % The stability of the quotient is typically determined by the stability
   % of the denominator. If stable, there is a reasonable expectation of
   % geometric decay in absolute value of the quotient coefficient
   % sequence. For well-conditioned denominators, the eigenvalue of the
   % largest magnitude (the "spectral radius") is often a good indication
   % of the number of lags required for the coefficient sequence to decay
   % to some small value before termination.
   
  [isDenominatorStable,eigenvalues] = isStable(B);

   if isDenominatorStable
       
      % A zero-degree denominator involves division by a constant, and the
      % degree of the quotient is that of the numerator (no 1000 degree
      % limit imposed):
      
      if isempty(eigenvalues)
         maxDegree = C.Degree;
         degree    = maxDegree;
      else
         maxDegree = min(ceil(log(relTol) / log(max(abs(eigenvalues)))), 1000);
         degree    = maxDegree + windowSize;
      end
   else
      maxDegree = min(ceil(log(realmax) / log(max(abs(eigenvalues)))), 1000);
      degree    = maxDegree;
   end

else
   degree = maxDegree + windowSize; % Allow for the termination window
end

% Initialize some parameters for performance and convenience, as well as a 
% convergence warning flag, but only for non-zero-degree denominators. If
% the degree of the denominator is zero, then just scale coefficients, and
% the termination algorithm is unnecessary.

p = B.Degree;                   % Degree of the denominator
isNonZeroDegree = p > 0;        % Is the denominator degree non-zero?
issueWarning = isNonZeroDegree; % Issue a convergence warning?

% Compute the quotient polynomial coefficients:

coefficients = cell(1,degree+1);     % Create an empty cell array
coefficients{1} = C.Coefficients{0} / B.Coefficients{0};
maximumNorm = norm(coefficients{1}); % Initialize largest norm so far
iCount = 0;                          % Number of consecutive lags <= termination tolerance

for i = 2:numel(coefficients)
    
    % The algorithm for polynomial division is derived from the convolution
    % formula designed to multiply two polynomials A(L) and B(L). See [2],
    % p. 373, equation 6.1.16, subject to the modifications discussed on
    % the bottom of p. 381 and the top of p. 382, and also allow for the
    % possibility that A0, B0, and C0 are not identity matrices.
    % 
    % To illustrate the approach, define A(L), B(L) and C(L) as polynomials
    % of degree p, d, and q, respectively:
    %
    %   A(L) = A(0) + A(1)*L + ... A(p)*L^p
    %   B(L) = B(0) + B(1)*L + ... B(p)*L^d
    %   C(L) = C(0) + C(1)*L + ... C(q)*L^q
    %
    % The modified equations are:
    %
    %   A(0)*B(0)                                                               = C(0)
    %   A(0)*B(1)   + A(1)*B(0)                                                 = C(1)
    %   A(0)*B(2)   + A(1)*B(1)   + A(2)*B(0)                                   = C(2)
    %       .                     .                        .                    =  .
    %       .                     .                        .                    =  .
    %       .                     .                        .                    =  .
    %   A(0)*B(p-1) + A(1)*B(p-2) + A(2)*B(p-3) + ... + A(p-1)*B(0)             = C(p-1)
    %   A(0)*B(p)   + A(1)*B(p-1) + A(2)*B(p-2) + ... + A(p-1)*B(1) + A(p)*B(0) = C(p)
    %   A(0)*B(p+1) + A(1)*B(p)   + A(2)*B(p-1) + ... + A(p-1)*B(2) + A(p)*B(1) = C(p+1)
    %   A(0)*B(p+2) + A(1)*B(p+2) + A(2)*B(p)   + ... + A(p-1)*B(3) + A(p)*B(2) = C(p+2)
    %       .                     .                        .                    =  .
    %       .                     .                        .                    =  .
    %       .                     .                        .                    =  .
    %
    % In the above sequence of equations, C(i) = 0 for all i > q.
    %
    % MRDIVIDE iteratively solves for the coefficients of A(L).

    coefficients{i} = C.Coefficients{i-1}; % Also handles 0-matrix initialization

    for j = 1:(i-1)
        if j <= p
           coefficients{i} = coefficients{i} - coefficients{i-j}*B.Coefficients{j};
        end
    end

    coefficients{i} = coefficients{i} / B.Coefficients{0};

%     if isNonZeroDegree % Nonzero degree denominator?
%         
%        % Assess termination criterion for nonzero degree denominators:
%        
%        currentNorm = norm(coefficients{i});
%        
%        % Apply the relative and absolute tolerance tests:
%        %
%        % o The norm of the current coefficient matrix is compared to the
%        %   largest coefficient norm found so far. If the ratio of the norms
%        %   is less than or equal to the relative tolerance, then the
%        %   relative termination criterion is satisfied.
%        %
%        % o If the magnitudes of all elements of the current coefficient
%        %   matrix are less than or equal to the absolute tolerance, then
%        %   the near-zero absolute termination criterion is satisfied.
%        
%        if ((currentNorm / maximumNorm) <= relTol) || all(abs(coefficients{i}(:)) <= absTol)
%            
%        % If this is the first lag of a window of consecutive lags that
%        % satisfies the termination tolerance, record the last lag whose
%        % coefficient was above tolerance. This identifies the beginning of
%        % a new tolerance window, and opens a window.
%        
%           if iCount == 0
%              degree = i-2;
%           end
%           
%           % Begin or continue counting the number of below-tolerance
%           % coefficients to satisfy the tolerance window:
%           
%           iCount = iCount+1;
% 
%           if iCount >= windowSize
%              issueWarning = false;
%              break % Tolerance window is complete, so stop
%           end
% 
%        else
% 
%           iCount = 0; % Reset the below-tolerance counter
% 
%           if currentNorm > maximumNorm
%              maximumNorm = currentNorm; % Record the new relative tolerance threshold
%           end
% 
%        end
% 
%     end
    
    % Check for nonfinite coefficients and numerical instabilities:
    
    if isnan(coefficients{i})
        
       warning('econ:LagOp:mrdivide:NaNEncountered',...
               'NaN-valued coefficients found. Terminating.')
           
       degree = i-2;
       issueWarning = false;
       break
    end

    if isinf(coefficients{i})
        
       warning('econ:LagOp:mrdivide:InfEncountered',...
               'Infinite-valued coefficients found. Terminating.')
           
       degree = i-2;
       issueWarning = false;
       break
    end

end

% Issue a warning about non-convergence, if necessary:

% if issueWarning
%    if iCount == 0
%        
%       wString = cat(2, 'Termination window not currently open and coefficients \n', ...
%                        '         are not below tolerance.');
%       warning('econ:LagOp:mrdivide:WindowNotOpen', wString)
%       
%    else
%       if iCount < windowSize
%           
%          warning('econ:LagOp:mrdivide:WindowIncomplete',...
%                  'Window component of termination \ncriterion not fully satisfied.')
%              
%       end
%    end
% end

warning(wState); % Restore warning state

% Enforce the maximum polynomial degree, apply the absolute tolerance, and 
% create the LagOp object:

degree = min(degree,maxDegree);
A = LagOp(coefficients(1:(degree+1)),'Lags',0:degree,'Tolerance',absTol);

end