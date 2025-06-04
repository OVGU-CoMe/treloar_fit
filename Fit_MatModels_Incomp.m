%% Fitting incompressible classic material models
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

% Define stress measurement (0 - 2nd Piola, 1 - 1st Piola, 2 - Cauchy)
flag.stress = 2;

% Define selected experiments; each row represents a fitting process
% [UT?, ET?, PS?]
%experim = [eye(3); 1,1,1];
experim = [1,1,1];

% Define filename
if strcmp(flag.matModel, 'Sw')
    filename = sprintf('Data_Fit_%s%s.mat', flag.matModel, ...
        num2str(flag.matNum));
else
    filename = sprintf('Data_Fit_%s.mat', flag.matModel);
end

%% Fitting process
% Set mandatory flag (only relevant for compressible models, not used here)
flag.J = 0;

% Loop through all data sets
for n1=1:size(experim,1)

    % Set flag for dataset
    flag.uniax = experim(n1,1);
    flag.equib = experim(n1,2);
    flag.pures = experim(n1,3);

    % Load selected dataset(s)
    if ~flag.uniax && ~flag.equib && ~flag.pures
        error('No experiment(s) selected!')
    else
        [lam1, lam2, lam3, S11] = deal([]);
    end

    if flag.uniax       % Load uniaxial tension data
        UT = load("Data_Treloar_UT.mat");
        lam1 = [lam1; UT.lam1];
        lam3 = [lam3; 1./sqrt(UT.lam1)];
        lam2 = [lam2; 1./sqrt(UT.lam1)];
        S11 = [S11; UT.S];
    end

    if flag.equib       % Load equibiaxial tension data
        ET = load("Data_Treloar_ET.mat");
        lam1 = [lam1; ET.lam1];
        lam3 = [lam3; ET.lam1];
        lam2 = [lam2; 1./ET.lam1.^2];
        S11 = [S11; ET.S];
    end

    if flag.pures       % Load pure shear data
        PS = load("Data_Treloar_PS.mat");
        lam1 = [lam1; PS.lam1];
        lam3 = [lam3; ones(size(PS.lam1))];
        lam2 = [lam2; 1./PS.lam1];
        S11 = [S11; PS.S];
    end

    % Load material model
    % Isochoric stress component in the 11-direction
    [matmod_iso_S11, ~, ~, numVar, lb, ub, con_A, con_b] = ...
        MaterialLaw(flag, 1, lam1, lam2, lam3, 3);
    % Isochoric stress component in the 22-direction
    [matmod_iso_S22, ~, ~, ~] = MaterialLaw(flag, 1, ...
        lam1, lam2, lam3, 4);
    % Volumetric stress component in the 11-direction
    matmod_vol_S11 = @(C) matmod_iso_S22(C).*(-1).*lam2.^(2-flag.stress).*lam1.^(flag.stress-2);
   
    %% Optimizer - ga with underlying fmincon
    tic
    % Options for local optimizer (fmincon)
    opt_local = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e+4);
    % Objective function for local optimizer (fmincon)
    Obj_local = @(C1) mse(S11.*lam1.^flag.stress, ...
        matmod_iso_S11(C1) + matmod_vol_S11(C1));
    % Optimized parameters using local optimizer
    C_opt_local = @(C) fmincon(Obj_local, C, con_A, con_b, ...
        [], [], lb, ub, [], opt_local);

    % Options for global optimizer
    opt_global = optimoptions('ga','PlotFcn', @gaplotbestf, ...
        'Generations', 50, 'PopulationSize', 20, 'UseParallel',true);
    % Objective function for global optimizer
    Obj_global = @(C) Obj_local(C_opt_local(C));
    % Genetic algorithm applied to objective function optimized by fmincon
    C_opt_global = ga(Obj_global, numVar, con_A, con_b,...
        [], [], lb, ub, [], [], opt_global);
    % Calculate Optimal parameters using fmincon function again
    C_opt = C_opt_local(C_opt_global);
   
    % Plot results
    figure;
    hold on
    plot(lam1, S11.*lam1.^flag.stress, 'DisplayName', 'Treloar Data')
    plot(lam1, matmod_iso_S11(C_opt)+matmod_vol_S11(C_opt), ...
        'DisplayName', 'Fitted model')
    legend('Location','northwest')

    %% Save material parameters
    if exist(filename, 'file') == 2
        Data = load(filename);
        Mat = Data.Mat;
    else
        Mat = struct;
    end

    if flag.uniax && ~flag.equib && ~flag.pures
        Mat.(sprintf('UT_%s', num2str(flag.stress))) = C_opt;
    elseif ~flag.uniax && flag.equib && ~flag.pures
        Mat.(sprintf('ET_%s', num2str(flag.stress))) = C_opt;
    elseif ~flag.uniax && ~flag.equib && flag.pures
        Mat.(sprintf('PS_%s', num2str(flag.stress))) = C_opt;
    elseif flag.uniax && flag.equib && flag.pures
        Mat.(sprintf('T_%s', num2str(flag.stress))) = C_opt;
    end
    Mat.UT_2=C_opt;
    Mat.ET_2=C_opt;
    Mat.PS_2=C_opt;
    Mat.T_2=C_opt;
    save(filename, "Mat")

end