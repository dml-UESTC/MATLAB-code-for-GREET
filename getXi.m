
function Xi = getXi(delta,Nt,Nr,L,s,thetaj,INR )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
pt = 1:Nt;
PT = -pt+pt';
pr = 1:Nr;
PR =  -pr+pr';
E = delta*log(kron(exp(PT),exp(PR)));
s = s/norm(s);
S = reshape(s,Nt,L);


Xi = eye(Nr*L);
for k = 1:length(thetaj)
    temp =   kron(S.',eye(Nr)) ;
    at = 1/sqrt(Nt)*exp(j*pi*(0:Nt-1)*(thetaj(k)))';
    ar = 1/sqrt(Nr)*exp(j*pi*(0:Nr-1)*(thetaj(k)))';
    Aj{k} = kron(at*at',ar*ar') ;
    Xi = Xi + 2/pi*10^(INR(k)/10)* temp * (  Aj{k} .*  sinc ( E) )* temp' ;
end





