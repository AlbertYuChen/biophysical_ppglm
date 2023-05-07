function [ktbasis,ktbasis_orig,ktbasprs] = make_Pillow_Basis(klength, nkbasis, b)
% ---- Set up fitting structure -------------------------------
gg.kt = [];       % basis weights for stimulus filter k
gg.ktbas = [];    % temporal basis for stimulus filter k
gg.ktbasis_orig = [];
gg.ktbasprs = []; % parameters for basis for k-filter

assert((klength>nkbasis), 'klength should be bigger than number of temporal basis vectors');

ktbasprs.neye = 0; % number of "identity" basis vectors
ktbasprs.ncos = nkbasis; % Number of raised-cosine vectors to use
ktbasprs.kpeaks = [0 klength*(1 - 1.5/nkbasis)];  % Position of 1st and last bump
ktbasprs.b = b; % Offset for nonlinear scaling (larger -> more linear), Originally it's 5
[ktbasis, ktbasis_orig] = makeBasis_StimKernel(ktbasprs,klength);
gg.ktbas = ktbasis;
gg.ktbasis_orig = ktbasis_orig;
gg.ktbasprs = ktbasprs;







