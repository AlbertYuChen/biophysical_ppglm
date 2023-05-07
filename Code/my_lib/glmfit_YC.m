function [bb,dev,stats,loglikeli,iter,b_stack] = glmfit_YC(x,y,distr,varargin)

if nargin < 2
    error(message('stats:glmfit:TooFewInputs'));
end

if nargin < 3 || isempty(distr)
    distr = 'normal';
else
    distr = lower(distr);
end

% Determine the syntax.
if nargin < 4
    newSyntax = false;
else
    arg = varargin{1};
    if ischar(arg) % either a link name (old syntax), or a parameter name
        try
            validatestring(arg, ...
                {'identity' 'log' 'logit' 'probit' 'comploglog' 'reciprocal' 'logloglink'});
            newSyntax = false;
        catch
            newSyntax = true;
        end
    else % power link exponent, or custom link, but not a parameter name
        newSyntax = false;
    end
end

% Process optional name/value pairs.
if newSyntax
    paramNames = {     'link' 'estdisp' 'offset' 'weights' 'constant' 'rankwarn' 'options' 'b0'};
    paramDflts = {'canonical'     'off'      []        []        'on'       true        [] []};
    [link,estdisp,offset,pwts,const,rankwarn,options,b0] = ...
                           internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

else % the old syntax glmfit(x,y,distr,link,estdisp,offset,pwts,const)
    link = 'canonical';
    estdisp = 'off';
    offset = [];
    pwts = [];
    const = 'on';
    rankwarn = true;
    options = [];
    b0 = [];
    if nargin > 3 && ~isempty(varargin{1}), link = varargin{1}; end
    if nargin > 4 && ~isempty(varargin{2}), estdisp = varargin{2}; end
    if nargin > 5 && ~isempty(varargin{3}), offset = varargin{3}; end
    if nargin > 6 && ~isempty(varargin{4}), pwts = varargin{4}; end
    if nargin > 7 && ~isempty(varargin{5}), const = varargin{5}; end
end

estdisp = internal.stats.parseOnOff(estdisp,'''estdisp''');

if isempty(options)
    iterLim = 100;
    convcrit = 1e-4;
else
    options = statset(statset('glmfit'),options);
    iterLim = options.MaxIter;
    convcrit = options.TolX;
end

% Separate y and N for binomial distribution
[y,N] = internal.stats.getGLMBinomialData(distr,y);

% Remove missing values from the data.  Also turns row vectors into columns.
[anybad,wasnan,y,x,offset,pwts,N] = statremovenan_YC(y,x,offset,pwts,N);
if anybad > 0
    switch anybad
    case 2
        error(message('stats:glmfit:InputSizeMismatchX'))
    case 3
        error(message('stats:glmfit:InputSizeMismatchOffset'))
    case 4
        error(message('stats:glmfit:InputSizeMismatchPWTS'))
%   case 5
        % N is empty, or was created from y (so its length must match)
    end
end

if isequal(const,'on')
    x = [ones(size(x,1),1) x];
end
dataClass = superiorfloat(x,y);
x = cast(x,dataClass);
y = cast(y,dataClass);

% Set distribution-specific information.
[link,estdisp,sqrtvarFun,devFun,muLims] = ...
    internal.stats.getGLMVariance(distr,estdisp,link,N,y,dataClass);

% If x is rank deficient (perhaps because it is overparameterized), we will
% warn and remove columns, and the corresponding coefficients and std. errs.
% will be forced to zero.
[n,ncolx] = size(x);
if isempty(pwts)
    [~,R,perm] = qr(x,0);
else
    [~,R,perm] = qr(x .* pwts(:,ones(1,ncolx)),0);
end
if isempty(R)
    rankx = 0;
else
    rankx = sum(abs(diag(R)) > abs(R(1))*max(n,ncolx)*eps(class(R)));
end
if rankx < ncolx
    if rankwarn
        warning(message('stats:glmfit:IllConditioned'));
    end
    perm = perm(1:rankx);
    x = x(:,perm);
else
    perm = 1:ncolx;
end

% Number of observations after removing missing data, number of coeffs after
% removing dependent cols and (possibly) adding a constant term.
[n,p] = size(x);

if isempty(pwts)
    pwts = 1;
elseif any(pwts == 0)
    % A zero weight means ignore the observation, so n is reduced by one.
    % Residuals will be computed, however.
    n = n - sum(pwts == 0);
end
if isempty(offset), offset = 0; end
if isempty(N), N = 1; end

% Instantiate functions for one of the canned links, or validate a
% user-defined link specification.
[linkFun,dlinkFun,ilinkFun] = stattestlink_YC(link,dataClass);

% Initialize mu
if isempty(b0)
    % Initialize mu and eta from y.
    mu = startingVals(distr,y,N);
    eta = linkFun(mu);
else
    % Initialize based on supplied coefficients
    if ~isvector(b0) || ~isreal(b0) || length(b0)~=size(x,2)
        error(message('stats:glmfit:BadIntialCoef'))
    end
    eta = offset + x * b0(:);
    mu = ilinkFun(eta);
end

% Set up for iterations
iter = 0;
warned = false;
seps = sqrt(eps);
b = zeros(p,1,dataClass);
loglikeli = 0;

b_stack = [];

while iter < iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    deta = dlinkFun(mu);
    z = eta + (y - mu) .* deta;
%     deta = 1./mu;
%     z = eta + y.* deta - ones(length(mu),1);

    % Compute IRLS weights the inverse of the variance function
%     sqrtirls = abs(deta) .* sqrtvarFun(mu);
%     sqrtw = sqrt(pwts) ./ sqrtirls;    
    sqrtirls = 1/sqrt(mu);
    sqrtw = sqrt(pwts.*mu);

    % If the weights have an enormous range, we won't be able to do IRLS very
    % well.  The prior weights may be bad, or the fitted mu's may have too
    % wide a range, which is probably because the data do as well, or because
    % the link function is trying to go outside the distribution's support.
    wtol = max(sqrtw)*eps(dataClass)^(2/3);
    t = (sqrtw < wtol);
    if any(t)
        t = t & (sqrtw ~= 0);
        if any(t)
            sqrtw(t) = wtol;
            if ~warned
                warning(message('stats:glmfit:BadScaling'));
            end
            warned = true;
        end
    end

    % Compute coefficient estimates for this iteration - the IRLS step
    b_old = b;
    loglikeli_old = loglikeli;
    eta_old = eta;
    mu_old = mu;
    
    [b,R] = wfit(z, x, sqrtw);

    % Form current linear predictor, including offset
    eta = x * b;

    % Compute predicted mean using inverse link function
    mu = ilinkFun(eta); 
%     mu = exp(eta); 

    % Force mean in bounds, in case the link function is a wacky choice
    if isscalar(muLims)
        % 'poisson' 'gamma' 'inverse gaussian'
        if any(mu < muLims(1))
            mu = max(mu,muLims(1));
        end
    elseif ~isempty(muLims)
        % 'binomial'
        if any(mu < muLims(1) | muLims(2) < mu)
            mu = max(min(mu,muLims(2)),muLims(1));
        end
    end
    
    btmp = zeros(ncolx,1,dataClass); btmp(perm) = b; b_stack = [b_stack btmp];
    
    % Check stopping conditions   
    loglikeli = y.*eta - mu;
    loglikeli = sum(loglikeli);
    converge_ind = abs(b) < max(abs(b))/20;
    
    disp(['IRT# ' num2str(iter) char(9)...
        'log-LL: ' num2str(loglikeli) char(9)...
        'D log-LL: ' num2str(abs(loglikeli-loglikeli_old))  char(9)... 
        'D |b_cov|_2: ' num2str(norm(b(converge_ind)-b_old(converge_ind))) char(9)... 
        'D |b_all|_2: ' num2str(norm(b-b_old)) char(9)... 
        '#Div: ' num2str(sum(converge_ind==0)) ...
        ])
    
    %---------- check for convergence of log likelihood ---------------
    if loglikeli < -1e10, break; end    
    if abs(loglikeli-loglikeli_old) < convcrit, break; end
    %---------- check for convergence of beta values ---------------
%     out_id = 24;
%     outliers = find(abs(b(out_id)-b_old(out_id)) > convcrit * max(seps, abs(b_old(out_id))) );
%     disp(['abs(b-b_old)'  num2str( abs(b(out_id)-b_old(out_id)) ) ])
%     disp(['threshold'  num2str( convcrit * max(seps, abs(b_old(out_id))) ) ])
    
%     if (~any(abs(b(converge_ind)-b_old(converge_ind)) > convcrit * max(seps, abs(b_old(converge_ind))) )), break; end 
%     if (~any(abs(b-b_old) > convcrit * max(seps, abs(b_old)) )), break; end
end
if iter > iterLim
    warning(message('stats:glmfit:IterationLimit'));
end

bb = zeros(ncolx,1,dataClass); bb(perm) = b; 
converge_ind(perm) = converge_ind;

disp('Diverge Index: ')
find(converge_ind' == 0)

if iter>iterLim && isequal(distr,'binomial')
    diagnoseSeparation(eta,y,N);
end

if nargout > 1
    % Sum components of deviance to get the total deviance.
    di = devFun(mu,y);
    dev = sum(pwts .* di);
end

% Return additional statistics if requested
if nargout > 2
    % Compute the sum of squares used to estimate dispersion, and the
    % Anscombe residuals.
    switch(distr)
    case 'normal'
        ssr = sum(pwts .* (y - mu).^2);
        anscresid = y - mu;
    case 'binomial'
        ssr = sum(pwts .* (y - mu).^2 ./ (mu .* (1 - mu) ./ N));
        t = 2/3;
        anscresid = beta(t,t) * ...
            (betainc(y,t,t)-betainc(mu,t,t)) ./ ((mu.*(1-mu)).^(1/6) ./ sqrt(N));
    case 'poisson'
        ssr = sum(pwts .* (y - mu).^2 ./ mu);
        anscresid = 1.5 * ((y.^(2/3) - mu.^(2/3)) ./ mu.^(1/6));
    case 'gamma'
        ssr = sum(pwts .* ((y - mu) ./ mu).^2);
        anscresid = 3 * (y.^(1/3) - mu.^(1/3)) ./ mu.^(1/3);
    case 'inverse gaussian'
        ssr = sum(pwts .* ((y - mu) ./ mu.^(3/2)).^2);
        anscresid = (log(y) - log(mu)) ./ mu;
    end

    % Compute residuals, using original count scale for binomial
    if (isequal(distr, 'binomial'))
        resid = (y - mu) .* N;
    else
        resid  = y - mu;
    end

    dfe = max(n - p, 0);
    stats.beta = bb;
    stats.dfe = dfe;
    if dfe > 0
        stats.sfit = sqrt(ssr / dfe);
    else
        stats.sfit = NaN;
    end
    if ~estdisp
        stats.s = 1;
        stats.estdisp = false;
    else
        stats.s = stats.sfit;
        stats.estdisp = true;
    end

    % Find coefficient standard errors and correlations
    if ~isnan(stats.s) % dfe > 0 or estdisp == 'off'
        RI = R\eye(p);
        C = RI * RI';
        if estdisp, C = C * stats.s^2; end
        se = sqrt(diag(C)); se = se(:);   % insure vector even if empty
        stats.covb = zeros(ncolx,ncolx,dataClass);
        stats.covb(perm,perm) = C;
        C = C ./ (se * se'); 
        stats.se = zeros(ncolx,1,dataClass); stats.se(perm) = se;
        stats.coeffcorr = zeros(ncolx,ncolx,dataClass);
        stats.coeffcorr(perm,perm) = C;
         
        lowerhalf = (tril(ones(p),-1)>0);
        rv = stats.coeffcorr(lowerhalf);
        Tstat = rv .* sqrt((dfe-2) ./ (1 - rv.^2));
        pp = zeros(p,class(x));
        pp(lowerhalf) = 2*tpvalue(-abs(Tstat),dfe-2);
%         stats.coeffcorr_pp = pp + pp' + diag(diag(stats.coeffcorr)); % Preserve NaNs on diag.
        stats.coeffcorr_pp = pp + pp' + eye(length(pp)); % Preserve NaNs on diag.
   
        stats.t = NaN(ncolx,1,dataClass); stats.t(perm) = b ./ se;
        if estdisp
            stats.p = 2 * tcdf(-abs(stats.t), dfe);
        else
            stats.p = 2 * normcdf(-abs(stats.t));
        end
    else
        stats.se = NaN(size(bb),class(bb));
        stats.coeffcorr = NaN(length(bb),class(bb));
        stats.t = NaN(size(bb),class(bb));
        stats.p = NaN(size(bb),class(bb));
        stats.covb = NaN(length(bb),class(bb));
    end

    stats.resid  = statinsertnan_YC(wasnan, resid);
    stats.residp = statinsertnan_YC(wasnan, (y - mu) ./ (sqrtvarFun(mu) + (y==mu)));
    stats.residd = statinsertnan_YC(wasnan, sign(y - mu) .* sqrt(max(0,di)));
    stats.resida = statinsertnan_YC(wasnan, anscresid);
    
    stats.wts = 1./sqrtirls.^2;
end


function [b,R] = wfit(y,x,sw)
% Perform a weighted least squares fit
[~,p] = size(x);
yw = y .* sw;
xw = x .* sw(:,ones(1,p));
% No pivoting, no basic solution.  We've removed dependent cols from x, and
% checked the weights, so xw should be full rank.
[Q,R] = qr(xw,0);
b = R \ (Q'*yw);


function mu = startingVals(distr,y,N)
% Find a starting value for the mean, avoiding boundary values
switch distr
case 'poisson'
    mu = y + 0.25;
case 'binomial'
    mu = (N .* y + 0.5) ./ (N + 1);
case {'gamma' 'inverse gaussian'}
    mu = max(y, eps(class(y))); % somewhat arbitrary
otherwise
    mu = y;
end


function diagnoseSeparation(eta,y,N)
% Compute sample proportions, sorted by increasing fitted value
[x,idx] = sort(eta);
if ~isscalar(N)
    N = N(idx);
end
p = y(idx);
if all(p==p(1))   % all sample proportions are the same
    return
end
if x(1)==x(end)   % all fitted probabilities are the same
    return
end

noFront = 0<p(1) && p(1)<1;     % no "front" section as defined below
noEnd = 0<p(end) && p(end)<1;   % no "end" section as defined below
if noFront && noEnd
    % No potential for perfect separation if neither end is perfect
    return
end

% There is at least one observation potentially taking probability 0 or
% 1 at one end or the other with the data sorted by eta. We want to see
% if the data, sorted by eta (x) value, have this form:
%
%            Front                Middle                End
%        ---------------     -----------------     ---------------
%        x(1)<=...<=x(A)  <  x(A+1)=...=x(B-1)  <  x(B)<=...<=x(n)
% with   p(1)=...=p(A)=0                           p(B)=...=p(n)=1
% or     p(1)=...=p(A)=1                           p(B)=...=p(n)=0
%        ---------------     -----------------     ---------------
%        x may vary here     x is constant here    x may vary here
%
% This includes the possibilities:
%     A+1=B  - no middle section
%     A=0    - no perfect fit at the front
%     B=n+1  - no perfect fit at the end
dx = 100*max(eps(x(1)),eps(x(end)));
n = length(p);
if noFront
    A = 0;
else
    A = find(p~=p(1),1,'first')-1;
    cutoff = x(A+1)-dx;
    A = sum(x(1:A)<cutoff);
end

if noEnd
    B = n+1;
else
    B = find(p~=p(end),1,'last')+1;
    cutoff = x(B-1)+dx;
    B = (n+1) - sum(x(B:end)>cutoff);
end

if A+1<B-1
    % There is a middle region with >1 point, see if x varies there
    if x(B-1)-x(A+1)>dx
        return
    end
end

% We have perfect separation that can be defined by some middle point
if A+1==B
    xmid = x(A) + 0.5*(x(B)-x(A));
else
    xmid = x(A+1);
    if isscalar(N)
        pmid = mean(p(A+1:B-1));
    else
        pmid = sum(p(A+1:B-1).*N(A+1:B-1)) / sum(N(A+1:B-1));
    end
end

% Create explanation part for the lower region, if any
if A>=1
    explanation = sprintf('\n   XB<%g: P=%g',xmid,p(1));
else
    explanation = '';
end

% Add explanation part for the middle region, if any
if A+1<B
    explanation = sprintf('%s\n   XB=%g: P=%g',explanation,xmid,pmid);
end
    
% Add explanation part for the upper region, if any
if B<=n
    explanation = sprintf('%s\n   XB>%g: P=%g',explanation,xmid,p(end));
end

warning(message('stats:glmfit:PerfectSeparation', explanation));

% ------------------------------------------------
function p = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic.

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
