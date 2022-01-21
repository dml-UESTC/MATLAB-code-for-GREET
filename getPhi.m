
function [ Phi ] = getPhi(delta,Nt,Nr,L,w,thetaj,INR )
pt = 1:Nt;
PT = pt-pt';
pr = 1:Nr;
PR =  pr-pr';
E = delta*log(kron(exp(PR),exp(PT))); 

W = reshape(w,Nr,L);
Phi = w'*w*eye(Nt*L);
for k = 1:length(thetaj)
    temp =   kron(W.',eye(Nt)) ;
    at = 1/sqrt(Nt)*exp(j*pi*(0:Nt-1)*(thetaj(k)))';
    ar = 1/sqrt(Nr)*exp(j*pi*(0:Nr-1)*(thetaj(k)))';
    Aj{k} = kron(conj(ar)*ar.',  conj(at)*at.') ;
    Phi = Phi + 2/pi*10^(INR(k)/10)* temp * (  Aj{k} .*  sinc ( E) )* temp' ;
end










