function t= update_t(r,s,rho1,u1,P,Gamr,PI)
       g = P'*(s-u1);
       eta = 1/(r'*Gamr*r);
       [nu]=bisection(g,PI,rho1,eta );
       gam = diag(PI);
       t_tilde =  rho1*g./( 2*eta*gam + rho1+2*nu );
       t = P*t_tilde;
end