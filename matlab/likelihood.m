% Edits by (tc): Modified theta to account for endmembers.
function L = likelihood(theta, D, LLK_prior)

% theta is a vector that has N grain sizes, N abundances, and N endmembers
% D is a structure containing the data D.rho, D.n, D.k, D_lam_SPEC, D.R_SPEC
theta = theta';

% MODEL ABUNDANCES AND GRAIN SIZES
theta = theta(theta ~= 0);
N = length(theta)/3; % # of endmembers
sizes = theta(1:N);
abundances = theta(N+1:2*N);
endmembers = theta(2*N+1:3*N);

% ENDMEMBERS SPECIFIC DATA
% Need to choose the right densities, wavelengths, and optical constants,
% depending on the endmembers selected by the algorithm.
densities = zeros(1, N);
lambda = {};
n = {};
k = {};
for i = 1:N
    index = endmembers(i);
    densities(1, i) = D.rho(index);
    lambda = [lambda D.lambda{index}];
    n = [n D.n{index}];
    k = [k D.k{index}];
end

% Get fractional surface area for mixing. Lapotre et al. (2017) assumes
% that the abundances are equivalent to the mass fraction (the first line).
% Note that other papers assume that the abundances are instead equivalent
% to the relative geometric cross section (the second line). This changes
% the answer, so need to be careful!
fi = wtperctof(densities, sizes, abundances);
%fi = abundances;

lam_SPEC = D.lam_SPEC;
R_SPEC = D.R_SPEC;
inc = D.inc;
e = D.e;
R_vs_SSA = D.R_vs_SSA;
lmin = D.lmin;
lmax = D.lmax;

% FORWARD MODEL
% Model the SSA of the mixture corresponding to theta
[lam_MIX, SSA_MIX, ~] = Hapke_opt_forward(endmembers, D.labels, fi, sizes, n, k, lambda);

% If the data spectrum is in units of reflectance, then convert the
% forward-modelled SSA spectrum to reflectance.
if R_vs_SSA
    mu0 = cos(inc*pi/180);      
    mu = cos(e*pi/180);
    R_MIX = reflect(SSA_MIX, mu0, mu);
else
    R_MIX = SSA_MIX;
end

% Calculate the likelihood.
ind = find(lam_SPEC >= lmin & lam_SPEC <= lmax);
R_int = interp1(lam_MIX, R_MIX, lam_SPEC(ind));
res = R_SPEC(ind) - R_int;
L = -0.5 * res' * D.Cinv * res;
