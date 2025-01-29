% Edits by (tc): Modified theta to account for endmembers.
function llk = llk_prior(theta, D)

% theta is a vector that has N grain sizes, N abundances, and N endmembers
% D is a structure containing the data D.rho, D.n, D.k, D_lam_SPEC, D.R_SPEC

% MODEL ABUNDANCES AND GRAIN SIZES
theta = theta(theta ~= 0);
N = length(theta)/3; % # of endmembers
sizes = theta(1:N);
abundances = theta(N+1:2*N);
endmembers = theta(2*N+1:3*N);

% Calculate the log-likelihood on basis of the grain size uniform PDF.
% (?) Why not on the basis of the abundance Dirichlet PDF as well? Also,
% should we consider the PDF of the randomly selected endmembers? If so,
% then how do we do this? Unlike the grain sizes, the endmembers are not
% drawn from a continuous probability distribution - it is discrete.
llk = 0;
for i = 1:N
    index = endmembers(i);
    llk = llk + log(unifpdf(theta(i), D.prior_size_low(index), D.prior_size_high(index)));
end

return
