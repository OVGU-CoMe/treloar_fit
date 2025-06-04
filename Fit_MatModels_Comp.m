%% Fitting compressible classic material models
%**************************************************************************
%   author: Lukas Maurer
%   mail:   lukas.maurer@ovgu.de
%   date:   28/01/2025
%
%**************************************************************************

clc
clear variables
close all

%% Define Fitting process

% Define Material model
% NH - Neo-Hooke, MR - Mooney-Rivlin, Is - Isihara, St - Steinmann,
% GT - Gent-Thomas, Sw - Swanson, Ye - Yeoh, AB - Arruda-Boyce,
% Ge - Gent, YF - Yeoh-Fleming, Ca - Carroll
flag.matModel = "Ye";

% Define number of material parameters (only for Swanson model)
flag.matNum = 8;

% Define initial bulk modulus K0
K0 = 500;

% Define stress measurement (0 - 2nd Piola, 1 - 1st Piola, 
%                               2 - mod. Cauchy, 3 - Cauchy)
flag.stress = 2;

% Define experiment type (UT - Uniaxial tension, ET - Equibiaxial tension
%                         PS - Pure shear, T - All experiments)
flag.exp = 'T';

%Set flag.J for stress calculation (Cauchy)
if flag.stress == 3
    flag.stress = 2;
    flag.J = 1;
else
    flag.J = 0;
end

% Define filenames for incompressible and compressible versions
if strcmp(flag.matModel, 'Sw')
    filename_com = sprintf('Data_Fit_%s%s_Comp_%s.mat', flag.matModel, ...
        num2str(flag.matNum), num2str(flag.stress+flag.J));
    filename_inc = sprintf('Data_Fit_%s%s.mat', flag.matModel, ...
        num2str(flag.matNum));
else
    filename_com = sprintf('Data_Fit_%s_Comp_%s.mat', flag.matModel, ...
        num2str(flag.stress+flag.J));
    filename_inc = sprintf('Data_Fit_%s.mat', flag.matModel);
end

% Define initial material parameters
if exist(filename_com,'file')
    load(filename_com)
    C0 = [C{end}(1:end-1), K0];
elseif exist(filename_inc,'file')
    load(filename_inc)
    C0 = [eval(sprintf('Mat.%s_%s', ...
        flag.exp, num2str(min(flag.stress,2)))), K0];
else
    error('No incompressible material parameters available!')
end

% Load dataset
lam1 = [];
S11 = [];
if strcmp(flag.exp,'UT') || strcmp(flag.exp,'T')    % Load UT
    UT = load("Data_Treloar_UT.mat");
    lam1 = [lam1; UT.lam1];
    S11 = [S11; UT.S];
end
if strcmp(flag.exp,'ET') || strcmp(flag.exp,'T')    % Load ET
    ET = load("Data_Treloar_ET.mat");
    lam1 = [lam1; ET.lam1];
    S11 = [S11; ET.S];
end
if strcmp(flag.exp,'PS') || strcmp(flag.exp,'T')    % Load PS
    PS = load("Data_Treloar_PS.mat");
    lam1 = [lam1; PS.lam1];
    S11 = [S11; PS.S];
end

% Load material model
% S for the isochoric part in 11 direction
[~, ~, ~, numVar, lb, ub, con_A, con_b] = ...
    MaterialLaw(flag, 1, lam1, lam1, lam1, 3);
lb = [lb, K0];
ub = [ub, K0];
if ~isempty(con_A)
    con_A = [con_A, zeros(size(con_A,1),1)];
end

tic
% Set options for local optimizer
opt_local = optimoptions('fmincon', 'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 1e+4);
% Define objective function for local optimizer
Obj_local = @(C) Optim_func_Compressible(C, lam1, flag, S11);
% Compute optimized parameters using local optimizer
[C_opt, resnorm] = fmincon(Obj_local, C0, con_A, con_b, [], [], lb, ub, ...
    [], opt_local);
toc

% Save material parameters from paper
if exist(filename_com, 'file') == 2
    load(filename_com);
    C = [C, C_opt];
    Delta = [Delta, resnorm];
else
    C = {C_opt};
    Delta = resnorm;
end
save(filename_com, 'C', 'Delta')