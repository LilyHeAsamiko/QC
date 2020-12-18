function H = Hamiltonian(p,q,r,s,d,theta, L, S1, S2)
    if d == 2
        H0 = -74.9749;
        XX = XXGate(theta, L);
        XX1 = XXGate(-theta, L);
        H = H0.*exp([-1,0;1,1].*conj(S1).*XX.*S1.*conj(S2).*XX1.*S2);
    end
end