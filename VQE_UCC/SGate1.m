function S1 = SGate1(L,r,d,theta)
    lambda = 1/r;
    S2 = 1;
    S1 = S2;
    p =2;
    q =2;
    r = 2;
    s = 2;    
    H = Hamiltonian(p,q,2,s,d,theta, L, S1, S2);
    for j = 1: L
        S1 = S1.*exp(-1i*theta*H.*lambda);
    end
end