function QSigma = QCOV(r)
    r = transpose(r);  %colum
    V = cov(r,r);
    l = 2;
    x = r(1:2:length(r)-1);%transpose(r(1:2:length(r)-1));
    p = r(2:2:length(r));%transpose(r(2:2:length(r)));
    B = [x,p];
    C = [ones(l,1), 1i*ones(l,1); ones(l,1), -1i*ones(l,1)];
    QSigma = 0.25*(B.*C)*V.*(conj(C).*conj(B))+0.5*ones(2*l,l); 
%    QSigma = 0.25*(B.*C)*V.*(conj(C).*conj(B))+0.5*ones(2*l,l); 
end