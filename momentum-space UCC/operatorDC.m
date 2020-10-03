function C_AB_IJ_qr = operatorDC(PX,PY,PZ,A,B,I,J,q,r)
    C_AB_IJ_qr = i/8;
    for P = B+1:A
        for Q = J+1:I%range(J+1,I)
            if r == 'Re'
                C_AB_IJ_qr = C_AB_IJ_qr.*PZ.*PZ.*(PY.*PY.*PY.*PX+PX.*PY.*PX.*PX+PY.*PY.*PX.*PY+PY.*PX.*PX.*PX-PY'.*PX.*PY.*PY-PX'.*PX.*PX.*PY-PX'.*PY.*PY.*PY-PX'.*PX.*PY.*PX);
            elseif r == 'Im'
                C_AB_IJ_qr = C_AB_IJ_qr.*PZ.*PZ.*(PY.*PY.*PY.*PY+PX.*PX.*PX.*PX+PX.*PY.*PX.*PY+PY.*PX.*PY.*PX+PX.*PY.*PY.*PX-PY'.*PX.*PX.*PY-PY'.*PY.*PX.*PX-PX'.*PX.*PY.*PY);
            end
        end
    end
end