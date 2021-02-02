%Instead of optimizing the subnetwork arange, we focus on getting the
%maximum cut of the two labels(manualy by two physicians) of EEGs and generate new labels. 
%pre-computed solution for intialization: The objective function is set as: 
%f(x1, x2, x3)=(N-p*q)^2 =128*x1*x2*x3-56*1*x2-48*1*x3+16*x2*x3-52*x1-52*x2-96*x3+196

%Reduce the 3-local term to 2-local term as follows: x1*x2*x3 = x4*x3 + 2*(x1*x2-2*x1*x4-2*x2*x4+3*x4) if x4 != x1*x2,:
%f(x1, x2, x3, x4) =200*x1*x2-48*x1*x3-512*x1*x4+16*x2*x3-512*x2*x4+128*x3*x4-52*x1-52*x2-96*x3+768*x4+196

%Use the positive of the coefficient of x1x2x3, to set the value for xi = (1-si)/2,  i= 1,2,3,4
%Get the result of f(x1, x2, x3, x4)=116*s1 + 100*s2 + 24*s3 -160*s4 +50*s1*s2 -12*s1*s3 -128*s1*s4 + 4*s2*s3 +4*s2*s3 -128*s2*s4 + 32*s3*s4 +298 = 2*g(s1,s2,s3,s4)

%Up date The energy function: HP(?z(i), ?z(j), ?z(i), ?z(j)) = g(s1,s2,s3,s4) = 58* ?z(1)+ 50*?z(2) + 12*?z(3) – 80* ?z(4) +25*?z(1)*?z(2)-6 *) + 80* ?z(4) +25*?z(1)?z(2)  -*6*?z(1) * ?z(3)-64* ?z(1) * ?z(4)+ 2*?z(2) * ?z(3)-64* ?z(2) * ?z(4) +16*?z(3) * ?z(4)+149*I

%The Ising function is to be optimized with the local fields: h’ =(?z(1) , ?z(2) , ?z(3) ,?z(4)) = (58, 50, 12, -80) 
%Couplings J = (?z(1) ,?z(2) , ?z(3) ,?z(4))’ , where ?z(1) = (0 ,25,-6,-64), ?z(2) = (0 ,0,2,-64), ?z(3) = (0,0,0,16), ?z(4) = (0,0,0,0) 

z = zeros(4,4);
z(1,:) = [0,25,-6,-64];
z(2,:) = [0,0,2,-64];
z(3,:) = [0,0,0,16];
z(4,:) = [0,0,0,0];
hi = [58;50;12;-80]; 
J = z';
%normalize z into {-1,+1}
zn = normalize(z,'range'); 

%Note that,aQUBO problem if xi?{0,1} and an Ising problem if xi?{?1,+1} 
T = 100; % fixed
t = 0:100;
ra = 1/100*t;
sA = ones(9, length(t));% for RA
QHGh = ones(9, 4);
QHGt = ones(9, 4, length(t));
QHGhpred = QHGh;
QHGtpred = QHGt;
for p = 0.1:0.9
    %RA
    sa = ones(size(t))*p;%*
    %decrease the anneal fraction stoapoint(tainv,sinv)
    for t = 0:40
        sa(t+1) = 1-0.4/200*t;
    end
    %allow for a pause,meaning we also allow a point (tbinv,sinv) at the same
    %sinv(*).Then increase aneal fraction again:
    for t = 80:100
        sa(t+1) = p+0.4/100*(t-400);
    end
    t = 0:100;
    sr = 10*(p-0.1):1/(length(t)-1):10*p;
    sA(10*p,:) = sa;
    %HG 
    %Hp here is the Q
    H0 = sum(hi.*z);
    Hp = zn;
    for j = 1:4
        for i = 1:max(j-1,1)
            Hp(i,j) = H0(i)+sum(J(i,j)*zn(i,:).*zn(j,:));
        end
    end
    %Qfinal with the number of points for random exploration is set to init points=100,th enumber of iterations for optimization is set to niter=200,and the noise level is set to alpha=0.01.
    % where no linear since alpha2 = 0 snd alpha2*z = 0
    alpha = 0.01;
    h = ones(length(hi),length(sr));
    for i = 1:4
            h(i,:) = sr*sum(-hi(i)*zn(i,:));
    end
    Qf = h;    
    for i = 1:length(sr)
        for j = 1:4
            for k = 1:max(j-1,1)
                Qf(j,:) = alpha*h(j,:)+Hp(k,j);
                QHGt(10*p,j,:) =  Qf(j,:);
            end
        end
    end
    QHGh(10*p,:) = sum(Qf, 2);
%     %thus
%     h = hi;
%     for i = 1:4
%         h(i,:) = sum(-hi(i)*zn(i,:));
%     end
end
figure()
plot(sA');
xlabel('time');
ylabel('s');
title('RA ');
figure()
plot(QHGh');
xlabel('hscale');
ylabel('p');
title('HG');

figure()
subplot(2,2,1)
pcolor(reshape(QHGt(:,1,:),size(QHGt,1),size(QHGt,3)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p and eigenvalues = ',num2str(QHGt(1,1,1))]);

subplot(2,2,2)
pcolor(reshape(QHGt(:,2,:),size(QHGt,1),size(QHGt,3)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(24)]);

subplot(2,2,3)
pcolor(reshape(QHGt(:,3,:),size(QHGt,1),size(QHGt,3)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(QHGt(1,3,1))]);

subplot(2,2,4)
pcolor(reshape(QHGt(:,4,:),size(QHGt,1),size(QHGt,3)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(QHGt(1,4,1))]);

figure()
subplot(2,2,1)
boxplot(reshape(QHGt(:,1,:),size(QHGt,3),size(QHGt,1)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p and eigenvalues = ',num2str(QHGt(1,1,1))]);

subplot(2,2,2)
boxplot(reshape(QHGt(:,2,:),size(QHGt,3),size(QHGt,1)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(24)]);

subplot(2,2,3)
boxplot(reshape(QHGt(:,3,:),size(QHGt,3),size(QHGt,1)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(QHGt(1,3,1))]);

subplot(2,2,4)
boxplot(reshape(QHGt(:,4,:),size(QHGt,3),size(QHGt,1)));
xlabel('t');
ylabel('p');
title(['HG time evolvment with p, eigenvalues = ',num2str(QHGt(1,4,1))]);
%Bayes + RA + HG 
tol = 0.1;
Res = ones(4,9,50);
for p = 0.1:0.9
    %RA
    sa = ones(size(t))*p;%*
    %decrease the anneal fraction stoapoint(tainv,sinv)
    for t = 0:40
        sa(t+1) = 1-0.4/200*t;
    end
    %allow for a pause,meaning we also allow a point (tbinv,sinv) at the same
    %sinv(*).Then increase aneal fraction again:
    for t = 80:100
        sa(t+1) = p+0.4/100*(t-400);
    end
    t = 0:100;
    sr = 10*(p-0.1):1/(length(t)-1):10*p;
    sA(10*p,:) = sa;
    %HG 
    %Hp here is the Q
    H0 = sum(hi.*z);
    Hp = zn;
    for j = 1:4
        for i = 1:max(j-1,1)
            Hp(i,j) = H0(i)+sum(J(i,j)*zn(i,:).*zn(j,:));
        end
    end
    %Qfinal with the number of points for random exploration is set to init points=100,th enumber of iterations for optimization is set to niter=200,and the noise level is set to alpha=0.01.
    % where no linear since alpha2 = 0 snd alpha2*z = 0
    alpha = 0.01;
    h = ones(length(hi),length(sr));
    for i = 1:4
            h(i,:) = sr*sum(-hi(i)*zn(i,:));
    end
    Qf = h;    
    QWf = h;  
    for i = 1:length(sr)
        W = ones(4,101);%max(j-1,1));
        for j = 2:4
            for k = 1:max(j-1,1)
                Qf(j,:) = alpha*h(j,:)+sum(J(j,k)*zn(j,:).*zn(k,:));
                res = sum(1000*abs(Qf(j,:)));
                %backward prediction
                for iters = 1:50 
                    QWftemp = alpha*h(j,:)+Hp(k,j) + sum(W(j,k)*zn(j,:).*zn(k,:)); 
                    temp = mean(zn(j,:).*zn(k,:)).*(1-2*Qf(j,:));
                    dW = sum(temp,1);
                    Wtemp = W(j,k)-dW; 
                    if mean(abs(Qf(j-1,:) - Wtemp*mean(zn(j,:).*zn(k,:))))< res
                        res = mean(abs(Qf(j,:) - W(j,:)*mean(zn(j,:).*zn(k,:))));
                        W(j,:) = Wtemp;
                        QWf(j-1,:) = QWftemp;
                    end
                    if res < tol
                        Res(j,10*p,iters) = res;
                        break
                    else
                        Res(j,10*p,iters) = Res(j,10*p,max(1,iters-1));
                    end
                end
            end
            QHGtpred(10*p,j,:) =  QWf(j-1,:);
        end
    end
    QHGhpred(10*p,:) = sum(QHGtpred(10*p,j,:), 3);
    %     %thus
%     h = hi;
%     for i = 1:4
%         h(i,:) = sum(-hi(i)*zn(i,:));
%     end
end
figure()
heatmap(W)
title('weight matrix')
figure()
plot(QWf')
title('origin HG')
HGpredict = mean(zn).*zn*W;
figure()
plot(100*HGpredict'-200)
title('predict HG')
Res2 = abs(QWf -HGpredict);
Res2norm = Res2./sum(Res2,1);
figure()
plot(Res2norm')
title('Residule')
KL = log(QWf/(HGpredict+0.00001)+0.00001)*QWf;

figure()
plot(KL)
title('KL converge')

zpredict = sqrt(abs(HGpredict/h));
meankl = kl;
stdkl = kl;
cv = kl;
for i = 2:length(kl)
    meankl(i) = mean(kl(1:i));
    stdkl(i) = std(kl(1:i));
    cv(i) = stdkl(i)/meankl(i);
end
figure()
subplot(2,1,1)
plot(cv)
title('CV')
subplot(2,1,2)
errorbar( meankl,stdkl)
title('KL errbar')
