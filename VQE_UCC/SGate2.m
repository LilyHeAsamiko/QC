function S2 = SGate2(L,r,d,theta)
    lambda = 1/r;
    S2 = 1;
    p =2;
    q =2;
    s = 2;    
    H = Hamiltonian(p,q,2,s,d,theta, L, S1, S2);
    for j = 1: L
        S2 = S2.*exp(-1i*theta*H*lambda/2);
    end
    for j = L:-1:1
        S2 = S2.*exp(-1i*theta*H*lambda/2);
    end
end