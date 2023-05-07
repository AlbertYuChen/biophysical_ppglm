
function F = Build_Link(nonLink_t, f_Link_NLN)

figure('Position', [200, 200, 1400, 500]);
subplot(121)
plot(nonLink_t, f_Link_NLN)
hold on
plot([-100 10], [-100 10], ':')


subplot(122)
plot(nonLink_t, f_Link_NLN)
hold on
plot([-100 10], [-100 10], ':')
ylim([-50 5])
xlim([-50 5])

% ----------- construct new link function ----
tRan = -20:0.01:5;
fRan = nonLink_t > -7 & nonLink_t < 2;
LF = @(x,xdata) x(1)*exp(x(2)*xdata + x(3)) + x(4);
x0 = [1 1 1 1];
[Coe,resnorm,~,exitflag,output] = lsqcurvefit(LF,x0, nonLink_t(fRan), f_Link_NLN(fRan));

% ----------- manually defined ------------
Coe = [-0.1453   -0.5648    1.1034   -0.1700];

% -----------------------
figure('Position', [200, 200, 600, 550]);
plot(nonLink_t, f_Link_NLN)
hold on
plot(tRan, LF(Coe, tRan))
ylim([-115 15])
xlim([-30 15])

figure('Position', [200, 200, 600, 550]);
plot(tRan, exp(LF(Coe, tRan)), 'Linewidth', 2 )
hold on
plot(tRan, exp(tRan) )
plot(tRan, 1 ./ (1+exp(-tRan)) )
plot(nonLink_t, 1 ./ (1+exp(-f_Link_NLN)) )
ylim([-1 5])
xlim([-6 3])


%% build link function
% link = @(mu) ( log( (log(mu)-Coe(4))/Coe(1) ) - Coe(3)) /Coe(2);
% dlink = @(mu) 1 ./ ( Coe(2) * mu .* (log(mu) - Coe(4)) );
% 
% dataClass = 'double';
% tiny = realmin(dataClass)^.25; % keep fourth powers from under/overflowing
% lowerBnd = log(tiny); upperBnd = -lowerBnd;
% ilink = @(eta) exp( constrain_YC(  Coe(1)*exp(Coe(2)*eta+Coe(3)) + Coe(4) ,lowerBnd,upperBnd) );

%% 
link = @(mu) linkExt(mu, Coe); 
dlink = @(mu) dlinkExt(mu, Coe);  
dataClass = 'double';
tiny = realmin(dataClass)^.25; % keep fourth powers from under/overflowing
lowerBnd = log(tiny); upperBnd = -lowerBnd;
ilink = @(eta) ilinkExt(constrain(eta,lowerBnd,upperBnd), Coe);

F = {link, dlink, ilink};
% save('F.mat', 'F')

t=[0.01:0.01:3]';
figure
plot(t, link(t))
% hold on
% plot(t, link0(t))
% plot(t, link1(t))
grid on

figure
plot(t, dlink(t))
% hold on
% plot(t, dlink0(t))
% plot(t, dlink1(t))
grid on

t = [-3:0.01:3]';
figure
plot(t, ilink(t))
% hold on
% plot(t, ilink0(t))
% plot(t, ilink1(t))
grid on





%% exp + straight line 
function out = linkExt(mu, Coe)
CC = -0.9;
out = ones(size(mu));
Ind1 = mu <= exp(CC);
Ind2 = mu > exp(CC);

% out = log(mu);
% out(Ind1) = (mu(Ind1) - exp(CC))/exp(CC) + CC;
% out = ( log( (log(mu)-Coe(4))/Coe(1) ) - Coe(3)) /Coe(2);
out(Ind1) = ( log( (log(mu(Ind1))-Coe(4))/Coe(1) ) - Coe(3)) /Coe(2);
out(Ind2) = log(mu(Ind2));


function out = dlinkExt(mu, Coe)
CC = -0.9;
out = ones(size(mu));
Ind1 = mu <= exp(CC);
Ind2 = mu > exp(CC);

% out(Ind1) = 1 / exp(CC);
out(Ind1) = 1 ./ ( Coe(2) * mu(Ind1) .* (log(mu(Ind1)) - Coe(4)) );
out(Ind2) = 1 ./ mu(Ind2);


function out = ilinkExt(eta, Coe)
CC = -0.9;
out = ones(size(eta));
Ind1 = eta <= CC;
Ind2 = eta > CC;

% out(Ind1) = (eta(Ind1)-CC)*exp(CC) + exp(CC);
out(Ind1) = exp(Coe(1)*exp(Coe(2)*eta(Ind1)+Coe(3)) + Coe(4));
out(Ind2) = exp(eta(Ind2));


function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;























