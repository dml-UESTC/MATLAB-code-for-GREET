function [mu]=bisection(g,PI,rho1,eta )

epsilon = 0.0000001;  
gam = diag(PI);

b=1000;
a=-eta*min(gam)-0.5*rho1 ;

k=1;
err=b-a;
while(abs(err)>epsilon)
    mu=(a+b)/2;
     f= func(mu, g,gam,rho1,eta );
    if(f>=0)
        a=mu;
    elseif(f<0)
        b=mu;
    end
    err=b-a;
    k=k+1;
end
