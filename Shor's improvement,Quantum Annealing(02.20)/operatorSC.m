function C_A_I_qr = operatorSC(PX,PY,PZ,A,I,q,r)
    C_A_I_qr = i/2;
    for M = q+1:A  %M = range(q+1,A)
        if r == 'Re'
            C_A_I_qr = C_A_I_qr.*PZ.*(PY'.*PX-PX.*PY);
        elseif r == 'Im'
            C_A_I_qr = C_A_I_qr.*PZ.*(PY'.*PX+PX.*PY);
        end
    end
end 
