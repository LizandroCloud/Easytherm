function [F] = pr(z,Tr,Pr,w)

    omega = 0.07780;
    beta = omega*Pr/Tr;
    psi = 0.45724;
    alpha = (1 + (0.37464+1.54226*w-0.26992*w^2)*(1-Tr^0.5) )^2;
    
    q = psi*alpha/(omega*Tr);
    e =  1-sqrt(2) ;
    lamb = 1+sqrt(2);
    F = z - (1 + beta - q*beta*(z-beta)/( (z+e*beta)*(z+lamb*beta) ));

