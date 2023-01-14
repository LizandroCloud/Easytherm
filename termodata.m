function [Tc,Pc,Vc,Zc,w] = termodata(molecule,state)
    Tc=0;
    %n-butane
    switch molecule
         case {'generic'}
            Tc = 190.6; %K
            Pc = 45.99e5; %Pa
            w = 0.012; %
            Vc = 220e-6; % m3
            Zc = 0.267; % dimensioneless
        case {'1-3-butadiene'}
            Tc = 425.2; %K
            Pc = 42.7*1e5; %Pa
            Vc = 220e-6; % m3
            Zc = 0.267; % dimensioneless
            w = 0.19    ; %
        case {'water'}
            Tc = 647.1; %K
            Pc = 220.55*1e5; %Pa
            Vc = 262e-6; % m3
            Zc = 0.282; % dimensioneless
            w = 0.345; %
        case {'i-butane'}
            Tc = 408.1; %K
            Pc = 36.48*1e5; %Pa
            Vc = 262e-6; % m3
            Zc = 0.282; % dimensioneless
            w = 0.181; %
         case {'methane'}
            Tc = 425.1; %K
            Pc = 37.96*1e5; %Pa
            Vc = 1; % m3
            Zc = 1; % dimensioneless
            w = 0; %
         case {'propane'}
            Tc = 369.8; %K
            Pc = 42.48*1e5; %Pa
            Vc = 1; % m3
            Zc = 1; % dimensioneless
            w = 0.152; %
         case {'ethylene'}
            Tc = 282.3; %K
            Pc = 50.4*1e5; %Pa
            Vc = 1; % m3
            Zc = 0.281; % dimensioneless
            w = 0.087; %
         case {'1-butene'}
            Tc = 420; %K
            Pc = 40.43*1e5; %Pa
            Vc = 262e-6; % m3
            Zc = 0.282; % dimensioneless
            w = 0.191; %
         case {'argon'}
            Tc = 150.9; %K
            Pc = 48.98*1e5; %Pa
            Vc = 74.6e-6; % m3
            Zc = 0.291; % dimensioneless
            w = 0; %
         case {'benzene'}
            Tc = 562.2; %K
            Pc = 48.98*1e5; %Pa
            Vc = 259e-6; % m3
            Zc = 0.271; % dimensioneless
            w = 0.21; %
         case {'n-butane'}
            Tc = 425.1; %K
            Pc = 37.96*1e5; %Pa
            Vc = 255e-6; % m3
            Zc = 0.274; % dimensioneless
            w = 0.2; %
         case {'carbon monoxyde'}
            Tc = 132.9; %K
            Pc = 34.99*1e5; %Pa
            Vc = 93.4e-6; % m3
            Zc = 0.299; % dimensioneless
            w = 0.048; %
         case {'carbon dioxyde'}
            Tc = 304; %K
            Pc = 73.83*1e5; %Pa
            Vc = 94e-6; % m3
            Zc = 0.274; % dimensioneless
            w = 0.224; %
        otherwise
    end
end 