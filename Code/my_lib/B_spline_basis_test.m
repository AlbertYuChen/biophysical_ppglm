% Author: Yu Chen
% Date: May 15th 2018 @ CNBC

% close all; clear; clc

%% 
% knots = [0, 0.01,  0.02, 0.03, 0.04, 0.05];
knots = [0 0.001 0.002 0.004 0.006 0.009 0.012  0.015 0.022 0.03 0.05]; 
rangeval = [0, 0.05];
norder = 4;
nbasis = norder + length(knots)-2; 
fdBasisParam2 = create_bspline_basis(rangeval,nbasis,norder,knots);

plot(fdBasisParam2)