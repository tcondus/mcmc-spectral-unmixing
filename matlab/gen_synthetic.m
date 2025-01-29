% Edits by (tc): Added component SSA spectra to the returned matrix. The
% fractional cross sections are also returned.
function [M, f] = gen_synthetic(AB, d, ind, R_vs_SSA, D)

if sum(AB) > 1, AB = AB./100; end

% Load optical constants for each laboratory endmember
% Need to choose the right densities, wavelengths, and optical constants,
% depending on the endmembers selected by the algorithm.
N = length(ind);
rho = zeros(1, N);
lambda = {};
n = {};
k = {};
for i = 1:N
    index = ind(i);
    rho(1, i) = D.rho(index);
    lambda = [lambda D.lambda{index}];
    n = [n D.n{index}];
    k = [k D.k{index}];
end
inc = D.inc;
e = D.e;

% Get fractional surface area for mixing. Lapotre et al. (2017) assumes
% that abundances are equivalent to the mass fraction (the first line).
% Note that other papers assume that the abundances are instead equivalent
% to the relative geometric cross section (the second line). This changes
% the answer, so need to be careful!
f = wtperctof(rho', d, AB);
%f = AB;

[lam, w, wif] = Hapke_opt_forward(ind, D.labels, f, d, n, k, lambda);

if R_vs_SSA
    mu0 = cos(inc*pi/180);
    mu = cos(e*pi/180);
    w = reflect(w, mu0, mu);
end

M = [lam' w wif];

breakpoint = 1;

end
