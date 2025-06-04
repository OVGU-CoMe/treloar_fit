%% Compute the fitting error of a hyperelastic material
%**************************************************************************
%   author: Lukas Maurer
%   mail:   lukas.maurer@ovgu.de
%   date:   28/01/2025
%
%**************************************************************************
%% Input
% flag:
%   .exp            Used Experiment
%                   UT - Uniaxial tension, ET - Equibiaxial tension
%                   PS - Pure shear, T - All experiments
%   .J:             Include J in the calculation of the stress?
%                   0 - No, 1 - Yes
%   .matModel:      Name of the hyperelastic material model
%                   NH - Neo-Hooke, MR - Mooney-Rivlin, Is - Isihara,
%                   St - Steinmann, GT - Gent-Thomas, Sw - Swanson,
%                   Ye - Yeoh, AB - Arruda-Boyce, Ge - Gent,
%                   YF - Yeoh-Fleming, Ca - Carroll
%   .matNum:        Number of material parameters (only for Swanson model) 
%   .stress:        Stress measure (0 - 2PK, 1 - 1PK, 2 - Cauchy)
%                   0 - No, 1 - Yes
% C:                Material parameters as vector
% lam1:             Principal stretch
% S11:              Stress from dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Error_fun:        Error function handle (wrt) the material parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Error_fun = Optim_func_Compressible(C, lam1, flag, S11)

% Calculate lambda2 & lambda3
[lam2, lam3] =  Sim_lam2(flag, C);

% Calculate Material response
[matmod_iso, matmod_vol, J, ~, ~, ~, ~, ~] = MaterialLaw(flag, C, ...
                                                lam1, lam2, lam3, 1);

% Calculate Error function
Error_fun = mse(S11.*lam1.^flag.stress.*J.^(-flag.J), ...
                matmod_iso + matmod_vol);

end