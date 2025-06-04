%% Compute principal stretches for compressible material formulations
%**************************************************************************
%   author: Lukas Maurer
%   mail:   lukas.maurer@ovgu.de
%   date:   28/01/2025
%
%**************************************************************************
%% Input
% flag:
%   .J:             Include J in the calculation of the stress?
%                   0 - No, 1 - Yes
%   .matModel:      Name of the hyperelastic material model
%                   NH - Neo-Hooke, MR - Mooney-Rivlin, Is - Isihara,
%                   St - Steinmann, GT - Gent-Thomas, Sw - Swanson,
%                   Ye - Yeoh, AB - Arruda-Boyce, Ge - Gent,
%                   YF - Yeoh-Fleming, Ca - Carroll
%   .matNum:        Number of material parameters (only for Swanson model) 
%   .stress:        Stress measure (0 - 2PK, 1 - 1PK, 2 - Cauchy)
% C:                Material parameters as vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% lam2, lam3:       Principal stretches
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lam2, lam3] =  Sim_lam2(flag, C)

% Initialize
lam2 = [];
lam3 = [];

% Load dataset
if strcmp(flag.exp,'UT') || strcmp(flag.exp,'T')    % Load UT
    UT = load("Data_Treloar_UT.mat");
end
if strcmp(flag.exp,'ET') || strcmp(flag.exp,'T')    % Load ET
    ET = load("Data_Treloar_ET.mat");
end
if strcmp(flag.exp,'PS') || strcmp(flag.exp,'T')    % Load PS
    PS = load("Data_Treloar_PS.mat");
end

% Initialize symbolic variables
syms lam1s lam2s lam3s

% Load stress formulation for material in 22 direction
[matmod_iso, matmod_vol, ~, ~, ~, ~, ~ , ~] = MaterialLaw(flag, C, ...
                                                lam1s, lam2s, lam3s, 2);

% Compute stress in 22 direction
S = matmod_iso + matmod_vol;

% Define termination criterion
eta_ter = 1e-5;

%% Uniaxial tension test
if strcmp(flag.exp,'UT') || strcmp(flag.exp,'T')
    % Initialize lam for experiment
    lam1_UT = UT.lam1;
    lam2_UT = ones(size(lam1_UT));
    % Loop over all loadsteps 
    for n1 = 2:length(lam1_UT)
        % Reset counter
        counter_it = 0;
        % Reset current criterion
        eta_cur = 1;
        % Use lam2 from last step
        lam2_UT(n1) = lam2_UT(n1-1);
        while eta_ter<eta_cur
            % Update timer
            counter_it = counter_it+1;
            % Numerical differentiation
            dS1 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_UT(n1), ...
                lam2_UT(n1)-1e-6, lam2_UT(n1)-1e-6}));
            dS2 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_UT(n1), ...
                lam2_UT(n1)+1e-6, lam2_UT(n1)+1e-6}));
            Sc = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_UT(n1), ...
                lam2_UT(n1), lam2_UT(n1)}));
            dS = (dS2-dS1)/2e-6;
            % Evaluate dlam2 (Sc -> 0)
            dlam2 = -Sc/dS;
            % Update current criterion
            eta_cur = abs(dlam2/lam2_UT(n1));
            % Update lam2
            lam2_UT(n1) = lam2_UT(n1)+dlam2;
            % Stop process if the model does not converge within 20 steps
            if counter_it>20
                eta_cur = 0;
                lam2_UT = 0*lam2_UT + 1;
            end          
        end
    end
    lam2 = [lam2; lam2_UT];
    lam3 = [lam3; lam2_UT];
end

%% Equibiaxial tension test
if strcmp(flag.exp,'ET') || strcmp(flag.exp,'T')
    % Initialize lam for experiment
    lam1_ET = ET.lam1;
    lam2_ET = ones(size(lam1_ET));
    % Loop over all loadsteps
    for n1 = 2:length(lam1_ET)
        % Reset counter
        counter_it = 0;
        % Reset current criterion
        eta_cur = 1;
        % Use lam2 from last step
        lam2_ET(n1) = lam2_ET(n1-1);
        while eta_ter<eta_cur
            % Update timer
            counter_it = counter_it+1;
            % Numerical differentiation
            dS1 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_ET(n1), ...
                lam2_ET(n1)-1e-6, lam1_ET(n1)}));
            dS2 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_ET(n1), ...
                lam2_ET(n1)+1e-6, lam1_ET(n1)}));
            Sc = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_ET(n1), ...
                lam2_ET(n1), lam1_ET(n1)}));
            dS = (dS2-dS1)/2e-6;
            % Evaluate dlam2 (Sc -> 0)
            dlam2 = -Sc/dS;
            % Update current criterion
            eta_cur = abs(dlam2/lam2_ET(n1));
            % Update lam2
            lam2_ET(n1) = lam2_ET(n1)+dlam2;
            % Stop process if the model does not converge within 20 steps
            if counter_it>20
                eta_cur = 0;
                lam2_ET = 0*lam2_ET + 1;
            end
        end
    end
    lam2 = [lam2; lam2_ET];
    lam3 = [lam3; lam1_ET];
end

%% Pure shear
if strcmp(flag.exp,'PS') || strcmp(flag.exp,'T')
    % Initialize lam for experiment
    lam1_PS = PS.lam1;
    lam2_PS = ones(size(lam1_PS));
    % Loop over all loadsteps
    for n1 = 2:length(lam1_PS)
        % Reset counter
        counter_it = 0;
        % Reset current criterion
        eta_cur = 1;
        % Use lam2 from last step
        lam2_PS(n1) = lam2_PS(n1-1);
        while eta_ter<eta_cur
            % Update timer
            counter_it = counter_it+1;
            % Numerical differentiation
            dS1 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_PS(n1), ...
                lam2_PS(n1)-1e-6, 1}));
            dS2 = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_PS(n1), ...
                lam2_PS(n1)+1e-6, 1}));
            Sc = double(subs(S, {lam1s, lam2s, lam3s}, {lam1_PS(n1), ...
                lam2_PS(n1), 1}));
            dS = (dS2-dS1)/2e-6;
            % Evaluate dlam2 (Sc -> 0)
            dlam2 = -Sc/dS;
            % Update current criterion
            eta_cur = abs(dlam2/lam2_PS(n1));
            % Update lam2
            lam2_PS(n1) = lam2_PS(n1)+dlam2;
            % Stop process if the model does not converge within 20 steps
            if counter_it>20
                eta_cur = 0;
                lam2_PS = 0*lam2_PS + 1;
            end
        end
    end
    lam2 = [lam2; lam2_PS];
    lam3 = [lam3; ones(size(lam1_PS))];
end