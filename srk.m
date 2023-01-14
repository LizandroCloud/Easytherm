function [F] = srk(z,Tr,Pr,w)

    omega = 0.08664;
    beta = omega*Pr/Tr;
    psi = 0.42748;
    alpha = (1 + (0.48+1.574*w-0.176*w^2)*(1-Tr^0.5) )^2;
    
    q = psi*alpha/(omega*Tr);
    e= 0 ;
    lamb = 1;
    F = z - (1 + beta - q*beta*(z-beta)/( (z+e*beta)*(z+lamb*beta) ));


end