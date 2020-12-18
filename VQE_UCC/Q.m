function Qfun = Q(r, Qsigma) % m clicks
    x = transpose(r(1:2:length(r)-1));
    p = transpose(r(2:2:length(r)));
    Alpha = x + p*1i;
    da = Alpha - mean(Alpha); 
    QSigma = QCOV(r);
    Qtemp = exp(-da*transpose(da)*inv(QSigma*transpose(QSigma))/2)/pi/abs(sqrt(det(QSigma*transpose(QSigma)))); 
    Qtemp(isnan(Qtemp)) = 0;
    Qtemp(isinf(Qtemp)) = 1;    
    Qfun = Qtemp;
end