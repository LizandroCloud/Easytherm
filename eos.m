
function [component] = eos(state,molecule,equation)

% Routine for Calculating Thermodynamic Properties
% Copyright (c) 2018 by UFF-TEQ-NEO - Universidade Federal Fluminense - Niteroi
% Revision: 1.0 $ $ Date: 2018/01/11 18:01 $
%          
% Lizandro de Sousa Santos (lizandrosousa@id.uff.br)
% or lizandrossantos@gmail.com

% INPUT::

% state: state.T -> actual temperature
%        state.P -> actual pressure
% molecule: molecule.name -> molecule name
% equation: state equation type -> (ig) Ideal Gas, 
%                                  (vdw)Van der Waals, 
%                                  (rk) Redlick Kwong,
%                                  (srk) Soave/Redlich/Kwong, 
%                                  (pr) Peng Robinson

[Tc,Pc,Vc,Zc,w] = termodata(molecule); % accessing the thermodynamic data library...


Tr = state.T*(Tc)^-1;  % Reduced Temperature
Pr = state.P*(Pc)^-1;  % Reduced Pressure
R = 8.31416; % Gas Constant
options = optimset('Display','iter', 'Tolfun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals', 1e4);

switch equation
    case 'vdw'
        % initial estimate vapour
        z0=1e-8; % compression factor of ideal gas
         [z,fval] = fsolve(@(z) vdw(z,Tr,Pr,w),z0,options);
        state.V(1) = z*R*state.T / (state.P);
        
        z0=(1/8)*Pr/Tr; % estimate for liquid (Smith Van Ness, 7th edition, page 72)
%         [z,fval] = fsolve(@(z) vdw(z,Tr,Pr,w),z0,options);
        state.V(2) = z*R*state.T / state.P;
    case 'rk'
        % initial estimate vapour
        z0=1; % compression factor of ideal gas
        [z,fval] = fsolve(@(z) rk(z,Tr,Pr),z0,options);
        state.V(1) = z*R*state.T / state.P;
        
        z0=(1/8)*Pr/Tr; % estimate for liquid (Smith Van Ness, 7th edition, page 72)
%         [z,fval] = fsolve(@(z) rk(z,Tr,Pr),z0,options);
        state.V(2) = z*R*state.T / state.P;
        
end


component.z = z;
component.V = state.V;


