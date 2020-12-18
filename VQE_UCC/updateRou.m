function rou = updateRou(V, r, l, labels, a) 
    hbar = 2;
    l =2;
    m = sum(sum(labels == 1));% m clicks
    N = 2^m; %steps between j-1 and j    
    %Va = ones(hbar*(l-1),hbar*(l-1));
    %Vb = ones(hbar,hbar);
    %Vab = ones(hbar*(l-1),hbar); 
    Va = V(1:hbar*(l-1),1:hbar*(l-1));
    Vb = V(1:hbar,1:hbar);
    Vab = V(1:hbar*(l-1),1:hbar); 
    ra = r(1:hbar*(l-1));
    rb = r(hbar*(l-1)+1:end);
    V = [Va, Vab; transpose(Vab), Vb];
    V = diag(V);
    rbn = reshape(rb(1:size(Vb(:))),[hbar,hbar]);
    Van = Va - Vab*inv(Vb + ones(hbar)).*rbn;
    Vn = [Van, Vab; transpose(Vab), Vb];
    rn = r(1:hbar*(l-1)) - Vab.*inv(Vb + ones(2,2)).*rbn;
    QSigma = QCOV(r);
    QSigman = QSigma*transpose(QSigma);
    O = ones(size(QSigman)) - inv(QSigman);
    P0 = [0, 0; 0, 1];
    P1 = [1, 0; 0, 1] - P0;
    T0 = (-1)^0/sqrt(det(ones(size(QSigman))-O))/sqrt(det(QSigman)); 
    T1 = (-1)^1/sqrt(det(ones(size(QSigman))-O))/sqrt(det(QSigman)); 
    q = hbar*exp(-transpose(rbn)*inv(Vb + ones(2,2))*reshape(r(1:size(Vb(:))),[2,2])./sqrt(det(Vb+ones(2,2))));
    q(isnan(q)) = 0.5;
    rou0 = round(abs(T0.*T1)); 
    rou = a;
    for i = 1:min(N,100)
        if isempty(a)
            ptemp = sum(diag(rou0*P0));
            if labels(i) == 0
                rou(i) = ptemp; 
            else
                rou(i) = sum(diag(rou0*ones(size(P1))*P1))/(1-ptemp); 
            end
        elseif length(a) == length(q) 
    %        p = sum(a.*q);
            On = ones(size(QSigman)) - inv(QSigman);
            T0 = (-1)^0/sqrt(det(ones(size(QSigman))-On))/sqrt(det(QSigman)); 
            T1 = (-1)^1/sqrt(det(ones(size(QSigman))-On))/sqrt(det(QSigman)); 
    %        qn = hbar*exp(-transpose(rbn)*inv(Vb + ones(2,2))*rn/sqrt(det(Vb+ones(2,2))));
            rou0n = round(abs(T0.*T1)); 
            if labels(i) == 1 
                rou(i) = sum(a(i,:)*rou0n(i)-q(i)*rou0n(i)); 
            else
                rou(i) = sum(a(i,:)*q(i)./p*rou0n(i));
            end
        end
    end
    rou(isnan(rou))= 0;
    rou(isinf(rou))= 1;
end