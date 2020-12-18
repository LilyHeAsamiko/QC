function XX = XXGate(theta, L)
    XX = exp(-1i*theta);    
    if L == 2
        XX = XX*exp(([0,1;1,0])^2/2);
    end
end