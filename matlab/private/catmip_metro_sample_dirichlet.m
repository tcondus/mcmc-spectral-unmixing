function [x, LLK, Naccept, Nreject] = catmip_metro_sample_dirichlet(x0, llk0, N0, Sm, beta, D, C)
% [theta_new,P,Naccept,Nreject] = catmip_metro_sample_dirichlet(theta,theta_lead,Nsamples,myfun,sigma);
% myfun returns log likelihood
% Edited so that second half of theta vector is explored using Dirichlet
% proposal PDF.
%
% Sarah Minson, May 27, 2014
% April 14, 2015: Bug fix to make proposal PDF symmetric.  Now q is a
% uniform Dirichlet distribution.
% Please cite:
% Minson, S. E., M. Simons, and J. L. Beck (2013), Bayesian inversion for finite fault earthquake source models I - theory and algorithm, Geophys. J. Int., 194(3), 1701-1726, doi:10.1093/gji/ggt180.
%
% Edits by (tc): Added support for endmembers.

Nparam_old = length(x0);
x0 = x0(x0 ~= 0);
Naccept = 0;
Nreject = 0;
Nparam = length(x0);

x = zeros(Nparam, N0);
x(:, 1) = x0;

LLK = zeros(3, N0);
LLK(:, 1) = llk0;

Sm0 = Sm;
Sm = eye(Nparam);
Sm(1:2, 1:2) = Sm0(1:2, 1:2);

z = mvnrnd(zeros(Nparam, 1), Sm, N0)';
z2 = dirichletrnd_exc(ones(Nparam/3), N0);
z3 = zeros(Nparam/3, N0);
[~, num_lib] = size(D.mineral_classes);
for i = 1:N0
    z3(:, i) = sort(randperm(num_lib, Nparam/3))';
end

for k = 1:N0-1 % sample posterior for this level of tempering
    y = x(:, k) + z(:, k);
    %y = round(y);
    y(Nparam/3+1:2*Nparam/3) = z2(:, k); % Replace middle third with Dirichlet sample
    y(2*Nparam/3+1:end) = z3(:, k); % Replace last third with (pseudo)-random endmember number

    px = LLK(1, k);
    [LLKy, error_code] = catmip_llk_post(y, beta, D, C);
    py = LLKy(1);

    r = exp(py-px);
    u = rand;

    if u <= r & ~error_code
        x(:, k+1) = y;
        LLK(:, k+1) = LLKy;
        Naccept = Naccept + 1;
    else
        x(:, k+1) = x(:, k);
        LLK(:, k+1) = LLK(:, k);
        Nreject = Nreject + 1;
    end
end

x_temp = zeros(Nparam_old, N0);
Nparam_diff = (Nparam_old / 3) - (Nparam / 3);
num_endmembers = Nparam / 3;
for i = 1:num_endmembers
    x_temp(i, :) = x(i, :); % x_temp is longer, x is shorter
    x_temp(num_endmembers+i+Nparam_diff, :) = x(num_endmembers+i, :);
    x_temp(2*num_endmembers+i+2*Nparam_diff, :) = x(2*num_endmembers+i, :);
end
x = x_temp;

return
