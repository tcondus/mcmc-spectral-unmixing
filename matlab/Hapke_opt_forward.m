% Edits by (tc): Allowed for grain sizes 10 microns in diameter or below to
% avoid a bug with smaller values. However, caution must be used because
% Hapke theory is only valid for the geometric optics regime (i.e., the
% grain size is "much larger" than the wavelength). The SSAs of the
% component spectra are also returned.
function [lam, w, wif] = Hapke_opt_forward(endmembers, labels, fi, Di, n, k, lambda)
%--------------------------------------------------------------------------
% CALCULATES REFLECTANCE OF A MIXTURE FROM SPECTRA OF INDIVIDUAL COMPONENTS
% Inputs:
%   1) fi = vector of relative cross-sections
%   2) Di = vector of grain sizes
% Dependences:
%   1) nkDtowi.m
%   2) reflect.m
%   3) Chandrasekhar.m
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

% Grain size in microns
%if Di(1)>10
if Di(1)>0
    Di = Di.*1e-6;
end

% Linear mixture of individual SSA
% Define wavelength over which the mixture spectrum will be modeled
% (the wavelength range is limited by the shortest range in the spectra)
lm = [];
lM = [];
N = length(fi);
for i = 1:N
    lm = [lm min(lambda{i})];
    lM = [lM max(lambda{i})];
end

lm = max(lm);
lM = min(lM);

dd = mean(diff(lambda{1}));
dl = dd/2;
lam = lm:dl:lM;

% Calculate SSA for each component of the mixture and interpolate it in the
% wanted wavelengths
for i = 1:N
    wi{i} = nkDtowi(n{i},k{i},Di(i),lambda{i});
	wif{i} = interp1(lambda{i},wi{i}',lam');
end

% Calculate SSA of the mixture
w = zeros(size(wif{1}));
for i = 1:N
    w = w+fi(i).*wif{i};
end
wif = cell2mat(wif); % Also return the SSAs of the component spectra.
