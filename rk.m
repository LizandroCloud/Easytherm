function [F] = rk(z,Tr,Pr,w)

    omega = 0.08664;
    beta = omega*Pr/Tr;
    psi = 0.42748;
    alpha = Tr^(-0.5);
    
    q = psi*alpha/(omega*Tr);
    e= 0 ;
    lamb = 1;
    F = z - (1 + beta - q*beta*(z-beta)/( (z+e*beta)*(z+lamb*beta) ));


end