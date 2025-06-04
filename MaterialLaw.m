%% Compute stress formulation for hyperelastic material law
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
% C1:               Material parameters as vector (only used if flagC<3)
% lam1, lam2, lam3: Principal stretches
% flagC:            Define output -> Stress direction & function handle
%                   1 - Stress in 11 direction (material parameters C1)
%                   2 - Stress in 22 direction (material parameters C1)
%                   3 - Stress in 11 direction, function handle @C
%                   4 - Stress in 22 direction, function handle @C
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% matmod_iso:	    Isochoric part of the stress function
% matmod_vol:       Volumetric part of the stress function
% J:                Calculated J values (Vector)
% numVar:           Number of material parameters
% lb:               Lower boundary for fitting process
% ub:               Upper boundary for fitting process
% con_A, con_B:     Linear inequality for fitting process A*x<=b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matmod_iso, matmod_vol, J, numVar, lb, ub, con_A, con_b] = ...
    MaterialLaw(flag, C1, lam1, lam2, lam3, flagC)

% Switch lam1 & lam2 if needed
if flagC == 2 || flagC == 4 
    [lam1, lam2] = deal(lam2, lam1);
end

% Invariants
I1 = lam1.^2 + lam2.^2 + lam3.^2;
I2 = lam1.^2.*lam2.^2 + lam2.^2.*lam3.^2 + lam1.^2.*lam3.^2;
I3 = lam1.^2.*lam2.^2.*lam3.^2;
J = sqrt(I3);

% Initalize nonl. constraint
con_A = [];
con_b = [];

%% Material models
if strcmp(flag.matModel, 'NH')      % Neo-Hooke
    % Isochoric part of the stress function
    matmod_iso = @(C) (C(1)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)))...
        .*lam1.^flag.stress.*J.^(-flag.J);
    % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 1;
    lb = -Inf;
    ub = -lb;
elseif strcmp(flag.matModel, 'MR')  % Mooney-Rivlin
    % Isochoric part of the stress function
    matmod_iso = @(C) (2*C(1)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)) + ...
        2*C(2)*J.^(-4/3).*(I1-lam1.^2-2/3*I2.*lam1.^(-2)))...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 2;
    lb = [0,0];
    ub = [Inf, Inf];
elseif strcmp(flag.matModel, 'Is')  % Isihara
    % Isochoric part of the stress function
    matmod_iso = @(C) (2*C(1)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)) + ...
        2*C(2)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)).*(J.^(-2/3).*I1-3)+...
        2*C(3)*J.^(-4/3).*(I1-lam1.^2-2/3*I2.*lam1.^(-2)))...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 3;
    lb = [-Inf, -Inf, -Inf];
    ub = -lb;
elseif strcmp(flag.matModel, 'St')  % Steinmann
    % Isochoric part of the stress function
    matmod_iso = @(C) (2*C(1)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)) + ...
        2*C(2)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)).*(J.^(-2/3).*I1-3).*...
        (J.^(-4/3).*I2-3)+ ...
        2*C(2)*J.^(-4/3).*(I1-lam1.^2-2/3*I2.*lam1.^(-2)).* ...
        (J.^(-2/3).*I1-3))...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 2;
    lb = [-Inf, -Inf];
    ub = -lb;
elseif strcmp(flag.matModel, 'GT')  % Gent-Thomas
    % Isochoric part of the stress function
    matmod_iso = @(C) (2*C(1)*J.^(-2/3).*(1-1/3*I1.*lam1.^(-2)) + ...
        2*C(2)./I2.*(I1-lam1.^2-2/3*I2.*lam1.^(-2))) ...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 2;
    lb = [-Inf, -Inf];
    ub = -lb;
    con_A = [-1, -1];
    con_b = 0;
elseif strcmp(flag.matModel, 'Sw')  % Swanson
    z1 = reshape(1:flag.matNum, [], 4);
    % Isochoric part of the stress function
    matmod_iso = @(C) 0;
    for n1 = 1:(flag.matNum/4)
        matmod_iso = @(C) matmod_iso(C) + ...
            sign(C(z1(n1, 1))).*(10.^(abs(C(z1(n1, 1)))-20)-(10^-20)).*...
            J.^(-2/3).* (J.^(-2/3).*I1/3).^C(z1(n1, 2)).* ...
            (1-1/3*I1.*lam1.^(-2)) + ...
            sign(C(z1(n1, 3))).*(10.^(abs(C(z1(n1, 3)))-20)-(10^-20)).*...
            J.^(-4/3).* (J.^(-4/3).*I2/3).^C(z1(n1, 4)) .* ...
            (I1-lam1.^2-2/3*I2.*lam1.^(-2));
    end
    matmod_iso = @(C) matmod_iso(C) .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = flag.matNum;
    if flag.matNum == 4
        lb = [-30, -20, -30, -20];
        ub = -lb;
    elseif flag.matNum == 8
        lb = [-30, -30, -20, -20, -30, -30, -20, -20];
        ub = -lb;
    end
elseif strcmp(flag.matModel, 'Ye')  % Yeoh
    % Isochoric part of the stress function
    matmod_iso = @(C) 2*J.^(-2/3).* ...
        (C(1) + 2*C(2)*(J.^(-2/3).*I1-3) + ...
        3*C(3)*(J.^(-2/3).*I1-3).^2).* ...
        (1-1/3*I1.*lam1.^(-2)) ...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 3;
    lb = [-Inf, -Inf, -Inf];
    ub = -lb;
elseif strcmp(flag.matModel, 'AB')  % Arruda-Boyce
    % Isochoric part of the stress function
    c = [1/2, 1/20, 11/1050, 19/7000, 519/673750];
    matmod_iso = @(C) 2*C(1)*J.^(-2/3).* ...
        (1:5).*(J.^(-2/3).*I1).^(0:4)*(c'./C(2).^(0:4)').* ...
        (1-1/3*I1.*lam1.^(-2)) ...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 2;
    lb = [0, 0];
    ub = [Inf, Inf];
elseif strcmp(flag.matModel, 'Ge')  % Gent
    % Isochoric part of the stress function
    matmod_iso = @(C) C(1)*J.^(-2/3).* ...
        (C(2)-3)./(C(2)-I1).* ...
        (1-1/3*I1.*lam1.^(-2)) ...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 2;
    lb = [0, 0];
    ub = [Inf, Inf];
elseif strcmp(flag.matModel, 'YF')  % Yeoh-Fleming
    % Isochoric part of the stress function
    matmod_iso = @(C) 2*J.^(-2/3).* ...
        (C(1).*exp(-C(2)*(J.^(-2/3).*I1-3)./(C(4)-3)) + ...
        C(3)*(C(4)-3)./(C(4)-J.^(-2/3).*I1)).* ...
        (1-1/3*I1.*lam1.^(-2)) ...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 4;
    lb = [0, 0, 0, 3];
    ub = [ Inf,  Inf,  Inf, Inf];
elseif strcmp(flag.matModel, 'Ca')  % Carroll
    % Isochoric part of the stress function
    matmod_iso = @(C) (2*J.^(-2/3).* ...
        (C(1)+4*C(2)*(J.^(-2/3).*I1).^3).* ...
        (1-1/3*I1.*lam1.^(-2)) + ...
        C(3)*J.^(-4/3)./(J.^(-4/3).*I2).^0.5.*...
        (I1-lam1.^2-2/3*I2.*lam1.^(-2)))...
        .*lam1.^flag.stress.*J.^(-flag.J);
     % Fitting parameter (number parameters, lower/upper boundary)
    numVar = 3;
    lb = [-Inf, -Inf];
    ub = -lb;
end

% Volumetric part of the stress function
matmod_vol = @(C) C(end) * (J.^2-J).*lam1.^(-2)...
    .*lam1.^flag.stress.*J.^(-flag.J);

% Insert material parameters in function handle
if flagC <= 2
    matmod_iso = matmod_iso(C1);
    matmod_vol = matmod_vol(C1);
end

end