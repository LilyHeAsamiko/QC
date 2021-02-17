%pre-computed solution for intialization: The objective function is set as: 
%f(x1, x2, x3)=(N-p*q)^2 =128*x1*x2*x3-56*1*x2-48*1*x3+16*x2*x3-52*x1-52*x2-96*x3+196
hG = 200*x(1)*x(2)-48*x(1)*x(3)-512*x(1)*x(4)+16*x(2)*x(3)-512*x(2)*x(4)+128*x(3)*x(4)-52*x(1)-52*x(2)-96*x(3)+768*x(4)+196;
%Reduce the 3-local term to 2-local term as follows: x1*x2*x3 = x4*x3 + 2*(x1*x2-2*x1*x4-2*x2*x4+3*x4) if x4 != x1*x2,:
%f(x1, x2, x3, x4) =200*x1*x2-48*x1*x3-512*x1*x4+16*x2*x3-512*x2*x4+128*x3*x4-52*x1-52*x2-96*x3+768*x4+196
hGqubo = 116*s(1) + 100*s(2) + 24*s(3) -160*s(4) +50*s(1)*s(2) -12*s(1)*s(3) -128*s(1)*s(4) + 4*s(2)*s(3) +4*s(2)*s(3) -128*s(2)*s(4) + 32*s(3)*s(4) +298;  
%Use the positive of the coefficient of x1x2x3, to set the value for xi = (1-si)/2,  i= 1,2,3,4
%Get the result of f(x1, x2, x3, x4)=116*s1 + 100*s2 + 24*s3 -160*s4 +50*s1*s2 -12*s1*s3 -128*s1*s4 + 4*s2*s3 +4*s2*s3 -128*s2*s4 + 32*s3*s4 +298 = 2*g(s1,s2,s3,s4)
%Up date The energy function: HP(?z(i), ?z(j), ?z(i), ?z(j)) = g(s1,s2,s3,s4) = 58* ?z(1)+ 50*?z(2) + 12*?z(3) – 80* ?z(4) +25*?z(1)*?z(2)-6 *) + 80* ?z(4) +25*?z(1)?z(2)  -*6*?z(1) * ?z(3)-64* ?z(1) * ?z(4)+ 2*?z(2) * ?z(3)-64* ?z(2) * ?z(4) +16*?z(3) * ?z(4)+149*I

%The Ising function is to be optimized with the local fields: h’ =(?z(1) , ?z(2) , ?z(3) ,?z(4)) = (58, 50, 12, -80) 
%Couplings J = (?z(1) ,?z(2) , ?z(3) ,?z(4))’ , where ?z(1) = (0 ,25,-6,-64), ?z(2) = (0 ,0,2,-64), ?z(3) = (0,0,0,16), ?z(4) = (0,0,0,0) 

%Aggregate 
%n = 15;
%k = 4;

%Instead of optimizing maximum cut of the two labels, we focus on getting the
%the subnetwork arange 

n = 4;
K = 3; 
z = zeros(1,4);
z(1) = 58;
z(2) = 50;
z(3) = 12;
z(4) = -80;
hi = z; 
J = zeros(4,4);
J(1,:) = hi;
J(2,:) = [0,25,-6,-64];
J(3,:) = [0,0,2,-64];
J(4,:) = [0,0,0,16];
%when n = 4
% assume fully connected, s in {-1?1}
sx = [1,1,1,1;1,1,1,-1;1,1,-1,1;1,-1,1,1;-1,1,1,1;1,1,-1,-1;1,-1,1,-1;-1,1,1,-1;1,-1,-1,1;-1,1,-1,1;-1,-1,1,1;1,-1,-1,-1;-1,1,-1,-1;-1,-1,1,-1;-1,-1,-1,1;-1,-1,-1,-1];
T = 1.57
kB = 1
h = 1
E = kB*T;
Hp = zeros(1,16);
for i= 1:16 
    Hp(i) = 58*sx(i,1)+50*sx(i,2)+12*sx(i,3)-80*sx(i,4)+25*sx(i,1)*sx(i,2)-6*sx(i,1)*sx(i,3)-64*sx(i,1)*sx(i,4)+2*sx(i,2)*sx(i,3)-64*sx(i,2)*sx(i,4)+16*sx(i,3)*sx(i,4)+149;
end

%s(t) insert 1us pause into 1us anneal time
ta = 100; % anealing time
tp = 100; %pause time: 100ns, 0.1us 
s = 50;
%ra = 1/100*t;
%example with pause:
t = [0:ta+tp-1];
sf = [0:s-1,s*ones(1,tp),s:ta-1]; 
sr = [ta-1:-1:ta-s,s*ones(1,tp),s:ta-1]; 
figure()
plot(sf/ta);
hold on 
plot(sr/ta,'--');
legend(['Forward';'Reverse'])
title('annealing parameter s as a functiion of time t')
xlabel('t(us)')
ylabel('s')
% without paulse
sfd = 0:ta-1;

% assume units T = 12.1mK = 1.57GHz, kB = 1, h = 1, Wcut = 1THz, eta = 10^(-3), L = 10
AB = xlsread('D:\Conference\Germany\09-1216A-A_DW_2000Q_6_annealing_schedule_data.xlsx');
sAB = AB(:,1);
A = AB(:,2);
B = AB(:,3);
C = AB(:,4);
ds = length(A)/length(sf);
dsd = length(A)/length(sfd);
%Iteratively correct HQ later or Alternatively, ransomize HQ 
%variational benchmark
[minHp,iid] = min(Hp);
H0 = -1/2*sum(sx(iid));
HQ = A(1:ds:end)*H0+B(1:ds:end).*repmat(Hp(iid),200,1);
HQd = A(1:dsd:end)*H0+B(1:dsd:end).*repmat(Hp(iid),100,1);

figure()
plot(A(1:ds:end),'b')
hold on 
plot(B(1:ds:end), 'r')
hold on
% plot([1 200],[E E], 'o')
% hline = gline;
% set(hline, 'k--');
yline(E, '--')
legend({'A(s)';'B(s)';'E = kBT'})
xlabel('s')
ylabel('Schedule(GHz)')
title('Annealing Schedule AB')

figure()
plot(0:1/199:1,A(1:ds:end)/B(1:ds:end)/200,'b')
hold on 
plot(C(1:ds:end)/200, 'r--')
legend({'Q(s)';'C(s) = kBT/B(s)'})
xlabel('s')
ylabel('Dimensionaless Schedule(GHz)')
title('Annealing Schedule QC')

%H = HQ+ HB +VQB 
%     = HT(x1,x2)(motion and trapping)+ Hi(t,x1,x2) = He(t,x1,x2)(HT-F) +
%       Hi(t,x1,x2)(avoid motional effects arising from He, through motion described by , He?t ? 1)
%     e splitting of the Hamiltonian according to Eq. (3b) to be meaningful we require the initial width of the atomic wave function a, as determined by
%     the trap, to be much smaller than the mean separation between the atoms R. We expand the dipole-dipole interaction around R, Vdip(ˆr) = u(R) ? F(ˆrz ? R) + . . .,
%     with F = 3u(R)/R. Here the first term gives the energy shift if both atoms are excited to state |ri, while
%     the second term contributes to He and describes the mechanical force on the atoms due to Vdip. Other contributions to He arise from the photon kick in the absorption |gi ? |ri,

%     non  = HQ +HLS-i/2*sum(Cad*Ca)

%Note that,aQUBO problem if xi?{0,1} and an Ising problem if xi?{?1,+1} 

Wc = 1*T;%THz
T = 1.57;%GHz
beta = 1/T;
eta = 10^(-3);
W = ones(4,3,200);%4*3*length(sf);
H = W;
g = W;
HB = sum(W,2);
phi = W;
VQB = W;
Norma = W;
Normad = W;
Normphi = W;
gamma = W;
Ca = W;
Nr = length(sf);%batch size
R = 100;
sigma = 1000;
Rp = 0;
resMC = min(HQ);
P = [];
% with paulse
while (resMC <= sigma) & (Rp <100)
    sigma = resMC;
    Rp = Rp +1;
    for i = 1:4
        for k = 1:3
            w = HQ*k*ds+ds;%H/k
            g(i,k,:) = sqrt(eta*w.*exp(-w/Wc)/(mean(w)*k-sum(w)));
            [Norma(:,k,:), Normad(:,k,:), Normphi(:,k,:)] = operatorA(sx(iid,:), k, ds, HQ);%phi(k+1)
            VQB(i,k,:) = g(i,k,:).*sum(z.^2*sum(Norma(i,k,:)+Normad(i,k,:)));
            gamma(i,k,:) = 2*pi*eta*w.*exp(-w/Wc)./(1-exp(-beta*w));
            Ca (i,k,:) = sqrt(gamma(i,k,:))*L;
            HB(i,:) = sum(w.*Normad(i,k,:).*Norma(i,k,:),2);
            H(i,k,:) = HQ(:,k)+ reshape(HB(i,:),200,1)+reshape(VQB(i,k,:),200,1);
            W(i,k,:) = w;
            p0 = ones(size(Normphi(i,k,:)))-Normphi(i,k,:).^2;
            resMC = min(abs(H(i,k,:)-HB(i,:)));
            Res(Rp) = resMC;
            if max(p0) < 0.5
                P(Rp)=1-(1-p0)^Rp;
%                 resMC = min(abs(H(i,k,:)-HB(i,:)));
%                 Res(Rp) = resMC;
            else
                break;
            end
        end
 end

% without paulse
Wd = ones(4,3,100);
Hd = Wd;
gd = Wd;
HBd = sum(Wd,2);
phid = Wd;
VQBd = Wd;
Normad = Wd;
Normadd = Wd;
Normphid = Wd;
gammad = Wd;
Cad = Wd;
Nrd = length(sfd);%batch size
Rd = 100;
sigmad = 1000;
Rpd = 0;
resMCd = min(HQ);
Pd = [];
while (resMCd <= sigmad) & (Rpd <100)
    sigmad = resMCd;
    Rpd = Rpd +1;
    for i = 1:4
        for k = 1:3
            wd = HQd*k*dsd+dsd;%H/k
            gd(i,k,:) = sqrt(eta*wd.*exp(-wd/Wc)/(mean(wd)*k-sum(wd)));
            [Normad(:,k,:), Normadd(:,k,:), Normphid(:,k,:)] = operatorA(sx(iid,:), k, dsd, HQd);%phi(k+1)
            VQBd(i,k,:) = gd(i,k,:).*sum(z.^2*sum(Normad(i,k,:)+Normadd(i,k,:)));
            gammad(i,k,:) = 2*pi*eta*wd.*exp(-wd/Wc)./(1-exp(-beta*wd));
            Cad (i,k,:) = sqrt(gammad(i,k,:))*L;
            HBd(i,:) = sum(wd.*Normadd(i,k,:).*Normad(i,k,:),2);
            Hd(i,k,:) = HQd(:,k)+ reshape(HBd(i,:),200,1)+reshape(VQBd(i,k,:),100,1);
            Wd(i,k,:) = wd;
            p0d = ones(size(Normphid(i,k,:)))-Normphid(i,k,:).^2;
            resMCd = min(abs(Hd(i,k,:)-HBd(i,:)));
            Resd(Rpd) = resMCd;
            if max(p0d) < 0.5
                Pd(Rpd)=1-(1-p0d)^Rpd;
            else
                break;
            end
        end
    end
end
      
%MCMC
% initialize with random s
Nr2 = length(sf);%batch size
R2 = 600;
sigma2 = 1000;
Rp2 = 1;
resMC2 = min(HQ);
P2 = zeros(1,R2);
W2 = ones(4,3,200);
H2 = W2;
g2 = W2;
HB2 = W2;
phi2 = W2;
VQB2 = W2;
Norma2 = W2;
Normad2 = W2;
Normphid2 = W2;
gammad2 = W2;
Ca2 = W2;

Nr2d = length(sfd);%batch size
R2d = 200;
sigma2d = 1000;
Rp2d = 1;
resMC2d = min(HQ);
P2d = ones(1,R2d);
W2d = ones(4,3,100);
H2d = W2d;
g2d = W2d;
HB2d = W2d;
phi2d = W2d;
VQB2d = W2d;
Norma1d = W2d;
Normad2d = W2d;
Normphid2d = W2d;
gammad2d = W2d;
Ca2d = W2d;

while Rp2 <= 600   
% with paulse
    sx02 = sx(randi(16),:);%1:16
    H02= -0.5*sum(sx02); 
    HQ2 = A(1:ds:end)*H02+B(1:ds:end).*repmat(Hp(iid),200,1);

    X2 = (1-sx(iid,:))/2;
%    Rp2 = Rp2 +1;   
    for i = 1:4
        for k = 1:3
            w2 = HQ2*k*ds+ds;%H/k
            g2(i,k,:) = sqrt(eta*w2.*exp(-w2/Wc)/(mean(w2)*k-sum(w2)));
            [Norma2(:,k,:), Normad2(:,k,:), Normphi2(:,k,:)] = operatorA(sx(iid,:), k, ds, HQ2);%phi(k+1)
            VQB2(i,k,:) = g2(i,k,:).*sum(z.^2*sum(Norma2(i,k,:)+Normad2(i,k,:)));
            gamma2(i,k,:) = 2*pi*eta*w2.*exp(-w2/Wc)./(1-exp(-beta*w2));
            Ca2 (i,k,:) = sqrt(gamma2(i,k,:))*L;
            HB2(i,k,:) = reshape(w2',1,200).*reshape(Normad2(i,k,:),1,200).*reshape(Norma2(i,k,:),1,200);
%            H2(i,k,:) = HQ2(:,k)+ reshape(sum(HB2(i,k,:),2),200,1)+reshape(VQB2(i,k,:),200,1);
            H2(i,k,:) = HQ2+ reshape(HB2(i,k,:),200,1)+reshape(VQB2(i,k,:),200,1);
            W2(i,k,:) = w2;
            p02 = ones(size(Normphi2(i,k,:)))-abs(Normphi2(i,k,:)).^2;
            mean(p02)
            Rp2
            if mean(p02) >= 0.5
                %P2d(Rp2d)=1-(1-p02d)^Rp2d;
                P2(Rp2)= mean(p02);%1-(1-max(p02))^Rp2;
                Rp2 = Rp2+1;
                %             else
%                 P2(Rp2:end)=[];
%                 break;
                resMC2 = mean(abs(H2(i,k,:)),2);
                if mean(abs(resMC2)) <= abs(sigma2)
                    Res2(Rp2) = mean(abs(resMC2))
                    sigma2 = mean(abs(resMC2));
                else
                    Res2(Rp2) = Res2(max(Rp2-1,1));
                end
            end
        end
    end
end
figure()
plot(HQ2)
title('HQ2')
HB2_re = reshape(sum(abs(HB2),2),size(HB2,1),size(HB2,3));
VQB2_re = reshape(sum(real(VQB2),2),size(VQB2,1),size(VQB2,3));
H2_re = reshape(sum(real(H2),2),size(H2,1),size(H2,3));

figure()
plot(H2_re')
xlabel('t')
ylabel('s')
title('MCMC H_k')
figure()
surf(HBd_re)
xlabel('t')
ylabel('s')
title('MCMC HB_k')
figure()
surf(VQB2_re')
xlabel('s')
ylabel('t')
title('MCMC VQB_k')
figure()
plot(P2(1:Rp2))
xlabel('Rp')
ylabel('P')
title('MCMC jump rate P')

PS= 1-(P2(1:Rp2).^Rp2);

while Rp2d < 200   
% without paulse
    HQ2d = A(1:dsd:end)*H02+B(1:dsd:end).*repmat(Hp(iid),100,1);
%    Rp2 = Rp2 +1;   
    for i = 1:4
        for k = 1:3
            w2d = HQ2d*k*dsd+dsd;%H/k
            g2d(i,k,:) = sqrt(eta*w2d.*exp(-w2d/Wc)/(mean(w2d)*k-sum(w2d)));
            [Norma2d(:,k,:), Normad2d(:,k,:), Normphi2d(:,k,:)] = operatorA(sx(iid,:), k, dsd, HQ2d);%phi(k+1)
            VQB2d(i,k,:) = g2d(i,k,:).*sum(z.^2*sum(Norma2d(i,k,:)+Normad2d(i,k,:)));
            gamma2d(i,k,:) = 2*pi*eta*w2d.*exp(-w2d/Wc)./(1-exp(-beta*w2d));
            Ca2d (i,k,:) = sqrt(gamma2d(i,k,:))*L;
            HB2d(i,k,:) = reshape(w2d',1,100).*reshape(Normad2d(i,k,:),1,100).*reshape(Norma2d(i,k,:),1,100);
%            H2d(i,k,:) = HQ2d(:,k)+ reshape(sum(HB2d(i,k,:),2),100,1)+reshape(VQB2d(i,k,:),100,1);
            H2d(i,k,:) = HQ2d+ reshape(HB2d(i,k,:),100,1)+reshape(VQB2d(i,k,:),100,1);
            W2d(i,k,:) = w2d;
            p02d = ones(size(Normphi2d(i,k,:)))-abs(Normphi2d(i,k,:)).^2;
            mean(p02d)
            Rp2d
            if mean(p02d) >= 0.5
                %P2d(Rp2d)=1-(1-p02d)^Rp2d;
                P2(Rp2d)= mean(p02d);%1-(1-max(p02d))^Rp2d;
                Rp2d = Rp2d+1;
                %             else
%                 P2(Rp2:end)=[];
%                 break;
                resMC2d = mean(abs(H2d(i,k,:)),2);
                if mean(abs(resMC2d)) <= abs(sigma2d)
                    Res2d(Rp2d) = mean(abs(resMC2d))
                    sigma2d = mean(abs(resMC2d));
                else
                    Res2d(Rp2d) = Res2d(max(Rp2d-1,1));
                end
            end
        end
    end
end   
HB2d_re = reshape(sum(abs(HB2d),2),size(HB2d,1),size(HB2d,3));
VQB2d_re = reshape(sum(real(VQB2d),2),size(VQB2d,1),size(VQB2d,3));
H2d_re = reshape(sum(real(H2d),2),size(H2d,1),size(H2d,3));

figure()
plot(H2d_re')
xlabel('t')
ylabel('s')
title('MCMC_noPaulse H_k')
figure()
surf(HB2d_re)
xlabel('t')
ylabel('s')
title('MCMC_noPausle HB_k')
figure()
surf(VQB2d_re')
xlabel('s')
ylabel('t')
title('MCMC_noPaulse VQB_k')
figure()
plot(P2d)
xlabel('Rpd')
ylabel('P2d')
title('MCMC noPausle jump rate')

figure()
boxplot(H2_re')
title('MCMC H Pulse')
figure()
boxplot(H2d_re')
title('MCMC H NoPulse')

figure()
boxplot(reshape(sum(real(H2),1),200,3))
title('MCMC H_k Pulse')
figure()
boxplot(reshape(sum(real(H2d),1),100,3))
title('MCMC H_k NoPulse')

figure()
bar([H2_re;H2d_re(1:2:end)]')
title('MCMC H')


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
