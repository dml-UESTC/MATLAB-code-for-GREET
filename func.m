function f= func(mu, g,gam,rho1,eta )

temp = rho1*g./( 2*eta*gam + rho1+2*mu);
f= temp'*temp-1;



