function s= update_s(r,t,rho1,rho2,u1,u2 )
    s =   1/(rho1+rho2)*( rho1*(t+u1)   +  rho2*(r+u2) ) ;
    am = 1/sqrt(length(s));
    for i=1:length(s)
        if    s(i) >= -am && s(i) <= am
           s(i)=s(i);
        elseif  s(i) < -am
            s(i)=-am;
        elseif s(i)>am
            s(i)=am;
        end
    end

