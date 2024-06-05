function [DR_pure,a1,a2,lambda1,lambda2,input_PDR,input_BAE,input_eta,mdl] = ...
    f_PDR_BAE_fitting(PDR_532,BSC_532,BSC_355,flag_dust,ixLayer,ixTime,BAE_pure)
% ----------------------------------------------------------------------- %
% SYLVA, https://sylva.bioaerosol.eu/
% @author: X. Shang
%
% Shang el al. 2022, https://acp.copernicus.org/articles/22/3931/2022/)
% non-linear least square regression fitting using lidar-derived 
% backscatter-related Ångström exponent and particle linear depolarization ratio
% ----------------------------------------------------------------------- %

if ~exist('BAE_pure','var')
    BAE_pure = 0;
end

BSC_355(:,flag_dust) = NaN;
BSC_532(:,flag_dust) = NaN;
PDR_532(:,flag_dust) = NaN;

BSC_355_layer = BSC_355(ixLayer,ixTime);
BSC_532_layer = BSC_532(ixLayer,ixTime);
PDR_532_layer = PDR_532(ixLayer,ixTime);

% ###### inputs for the fitting
lambda1 = 532;
lambda2 = 355;

% # Reshape the matrix into a 1D array
bsc1 = reshape(BSC_532_layer,1,[]);
bsc2 = reshape(BSC_355_layer,1,[]);
pdr1 = reshape(PDR_532_layer,1,[]);

% # Calculate BAE (backscatter related Angstrom exponent)
bae = - log(bsc1./bsc2) / log(lambda1/lambda2);
% # Calculate eta from bae
eta = (lambda2/lambda1).^(-bae);

% ###### Generate input values (remove NaN, remove Inf)
% # Create a boolean mask for NaN and Inf values
nan_mask = isnan(eta.*pdr1);
inf_mask = isinf(eta.*pdr1);

% # Combine the masks to identify NaN and Inf values
invalid_mask = nan_mask | inf_mask;

% # Use the mask to filter the array and remove NaN and Inf values
input_PDR = pdr1(~invalid_mask);
input_eta = eta(~invalid_mask);
input_BAE = bae(~invalid_mask);

% ###### Define a non-linear function to fit
% # Perform non-linear fitting
modelfun = @(a,x) (a(1)+a(2))./(x+1) - a(2);
mdl = fitnlm(input_PDR,input_eta,modelfun,[3 5]);

a1 = mdl.Coefficients.Estimate(1);
a2 = mdl.Coefficients.Estimate(2);

%% ###### calculate the DR (Depolarization Ratio) of pure pollen
DR_pure  = (a1+a2)/((lambda2/lambda1).^(-BAE_pure)  + a2) -1;


