% Author: Yu Chen
% Date: Aug 14th 2018 @ Beijing

function LL = log_likelihood(X, y, beta, distr, varargin)

if nargin > 2
    distr = convertStringsToChars(distr);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error('Too few inputs');
end

if nargin < 3 || isempty(distr)
    distr = 'bernoulli';
else
    distr = lower(distr);
end

% Process optional name/value pairs.
paramNames = {'offset'};
paramDflts = {      []};

option_args=[paramNames; paramDflts];
options = struct(option_args{:});

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[])
   inpName = lower(pair{1});

   if any(strcmp(inpName,paramNames)) 
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end


if isempty(options.offset)
    offset = 0;
else
    offset = options.offset;
end

N = length(y);

%% log-likelihood
eta = X*beta+offset;
second_term = 1 + exp(eta);
second_term = log(second_term);

log_likelihood = y' * eta - sum( second_term );
LL = log_likelihood;


























