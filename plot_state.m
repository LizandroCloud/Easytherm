
% Routine for Plotting a PVT State Equation
% Copyright (c) 2018 by UFF-TEQ-NEO - Universidade Federal Fluminense - Niteroi
% Revision: 1.0 $ $ Date: 2018/01/11 19:35 $
%          
% Lizandro de Sousa Santos (lizandrosousa@id.uff.br)
% By NEO - Núcleo de Estudos de Otimização
% https://neouff.wixsite.com/home

% INPUT

% state: state.T -> actual temperature
%        state.P -> actual pressure
% molecule: molecule.name -> molecule name
% equation: state equation type -> (ig) Ideal Gas, 
%                                  (vdw)Van der Waals, 
%                                  (rk) Redlick Kwong,
%                                  (srk) Soave/Redlich/Kwong, 
%                                  (pr) Peng Robinson

function [] = plot_state(tspan,pspan,molecule,equation)

    R = 8.31416; % Universal Gas Constant 
    N = tspan(end); % T range
    M = pspan(end); % P range
    rangeT = linspace(tspan(1),tspan(2),N);
    rangeP = linspace(pspan(1),pspan(2),M);
   
    switch equation
        case 'vdw'
        % Estimating the minimum volume
         state.T = rangeT(1); % actual T
         state.P = rangeP(M); % actual P
         [component] = eos(state,molecule,equation);  % computing v
         v_min = component.V(2); % for liquid
         
         % Estimating the maximum volume
         state.T = rangeT(N); % actual T
         state.P = rangeP(1); % actual P
         [component] = eos(state,molecule,equation);  % computing v
         v_max = component.V(1); % for gas
    end
    
     rangev = linspace(v_min,v_max*5,M);
     [Tc,Pc,Vc,Zc,w] = termodata(molecule); 
     
    for i=1:N
           for j=1:M 
                state.T = rangeT(i); % actual T
                state.v = rangev(j); % actual v
                    b = R*Tc / (8*Pc);
                    a = 27*R^2*Tc^2 / (64*(Pc));
                    state.P = (R*state.T / (state.v - b)) - (a / (state.v)^2);
               
%                 [component] = eos(state,molecule,equation);  % computing v
                
                v(i,j) = rangev(j);
                T(i,j) = rangeT(i);
                P(i,j) = state.P;
           end
           
           plot(v(i,:),P(i,:));
           hold on;
    end     

    figure(2);
    surf(v,T,P);
    
end