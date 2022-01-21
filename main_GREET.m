clc
clear
% Description: GREET for joint one-bit transmit waveform and receive filter design
% MATLAB version: 2016b

%-------------------parameter setting-------------------%
Nt = 10;
Nr = 8;
L =16;
theta0 = sind(34);
thetaj = sind([-50 10]);        % normalized angle

INR = 30*ones(size(thetaj));
SNR = [20];
delta = 0;        % uncertainty of inteference

rho1 = 2;
rho2 = 30;
num_AltOpt = 50;
num_ADMM = 200;

at0 = 1/sqrt(Nt)*exp(j*pi*(0:Nt-1)*(theta0))';
ar0 = 1/sqrt(Nr)*exp(j*pi*(0:Nr-1)*(theta0))';
A0 = kron(eye(L),ar0*at0.');

s = sign(randn(Nt*L,1))+ j* sign (randn(Nt*L,1)); 
s = s/norm(s);
t = 1/sqrt(2*Nt*L)* sign(randn(2*Nt*L,1));
r = 1/sqrt(2*Nt*L)* sign(randn(2*Nt*L,1));      % initialize variables for GREET

tic
num = 0;
while(num<num_AltOpt )
    num=num +1;
    %-------------------optimize w-------------------%
    Xi =  getXi(delta,Nt,Nr,L,s,thetaj,INR );
    Xi_inv = inv(Xi);
    w = Xi_inv *A0*s/(s'*A0'*Xi_inv *A0*s); % MVDR output
    
    %-------------------optimize s-------------------%
    Gam = A0'*w*w'*A0;
    Phi = getPhi(delta,Nt,Nr,L,w,thetaj,INR );
    Gamr = [ real(Gam)  -imag(Gam);  imag(Gam)  real(Gam)   ];
    Phir = [ real(Phi)  -imag(Phi);  imag(Phi)  real(Phi)   ];   % transform into real-valued form
    
    [U,Lam]=eig((Gamr+Gamr.')/2);
    [P,PI]=eig((Phir+Phir.')/2);            % EVD
    
    s_1bit = sign(real(s)) + j* sign(imag(s));   
    s_1bit = s_1bit/norm(s_1bit);               % ensure the objective recorded is computed by one-bit transmit waveform
    SINR(num) =  10*log10(  10^(SNR/10)  *  s_1bit'*Gam*s_1bit /(s_1bit'*Phi*s_1bit ));     % record objective

    u1 = zeros(2*Nt*L,1);
    u2 = zeros(2*Nt*L,1);        % initialize dual variables for ADMM
   
   for ii = 1:num_ADMM
        s = update_s(r,t,rho1,rho2,u1,u2 );               
        t = update_t(r,s,rho1,u1,P,Gamr,PI);
        r = update_r(t,s,u2,rho2,U,Phir,Lam);
        u1 = u1 + (t-s);
        u2 = u2 + (r-s);       
        
       er1(num,ii) = norm(s-t)^2;
       er2(num,ii) = norm(s-r)^2;
       er3(num,ii) =  abs(s'*s-1);
       er4(num,ii) =  norm(abs(s) - 1/sqrt(Nt*L*2) ); % record errors 
   end

    s = s(1:Nt*L)+j*s(Nt*L+(1:Nt*L)); % transform s into real-valued form
    fprintf('iter_num£º%d£¬SINR£º%5.4f dB£¬QSINR£º%5.4f dB£¬ ||s-t||= %5.4d£¬||s-r||= %5.4d  \n',num,SINR(num),SINR(num)+10*log10(2/pi),er1(end),er2(end) ) 
end
toc

figure(1)
QSINR = abs(SINR+10*log10(2/pi));
plot(1:length(SINR),QSINR,'-','LineWidth',2)
hold on
grid on
xlabel('outer-loop iteration')
ylabel('QSINR (dB)')

