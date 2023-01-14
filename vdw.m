function [F] = vdw(z,Tr,Pr,w)

    omega = 1/8;
    beta = omega*Pr/Tr;
    psi = 27/64;
    alpha = 1;
    q = psi*alpha/(omega*Tr);
    e= 0 ;
    lamb = 0;
    F = z - ( 1 + beta - q*beta*(z-beta)/( (z+e*beta)*(z+lamb*beta) ) );


end