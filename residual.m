% Subroutine for residual properties

function [Hr,Sr] = residual(state,molecule,equation)

% Routine for Calculating Thermodynamic Properties (Residual Properties)
% Copyright (c) 2018 by UFF-TEQ-NEO - Universidade Federal Fluminense - Niteroi
% Revision: 1.0 $ $ Date: 2018/01/11 17:35 $
%          
% Lizandro de Sousa Santos (lizandrosousa@id.uff.br)
% By NEO - N�cleo de Estudos de Otimiza��o
% https://neouff.wixsite.com/home

% INPUT

% state: state.T -> actual temperature
%        state.P -> actual pressure
% molecule: molecule.name -> molecule name (see termodata)
% equation: state equation type -> ig for Ideal Gas, 
%                                  vdw for Van der Waals, 
%                                  rk for Redlick Kwong,
%                                  srk for Soave/Redlich/Kwong, 
%                                  pr for Peng Robinson
%                                  gc for generalized correlations   



[Tc,Pc,Vc,Zc,w] = termodata(molecule,state); % accessing the thermodynamic data library...

Tr = (state.T)*(Tc)^-1;  % Reduced Temperature
Pr = (state.P)*1e5*(Pc)^-1;  % Reduced Pressure
R = 8.31416; % Gas Constant

switch equation
    
     case 'gc'
        B0  = 0.083 - 0.422/( (Tr)^1.6 );
        B1  = 0.139 - 0.172/( (Tr)^4.2 );
        dB0 = 0.675/( (Tr)^2.6 );
        dB1 = 0.722/( (Tr)^5.2 );

        Hr = R*Tc*Pr*( (B0 - Tr*dB0 + w*(B1 - Tr*dB1 )) );
        Sr = - R*Pr* (dB0 + w*(dB1) );
    case 'pr'
        P=state.P*1e5;
        T=state.T;
        om = 0.37464+1.54226*w-0.26992*w^2;
        b = 0.07780*R*Tc/Pc;
        a = (0.45724*R^2*Tc^2/Pc)*(1+om*(1-Tr^0.5))^2;
        af = a*P/(R*T)^2;
        beta=(b*P/(R*T));
        p = [1 (beta-1) af-3*beta^2-2*beta -af*beta+beta^3+beta^2];
        r = roots(p);
        z=min(r);
        da_dt = -0.45724*(R^2*Tc/Pc)*( (1+om*(1-sqrt(Tr)))*om )/(sqrt(Tr));
        Hr = R*T*(z-1) + ((T*(da_dt)-a) / (2*sqrt(2)*b))*log( ( ( 1+sqrt(2) )*beta + z )/((1-sqrt(2))*beta + z));
        Sr=0;
end
end
