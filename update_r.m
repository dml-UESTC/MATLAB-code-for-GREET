function r= update_r(t,s,u2,rho2,U,Phir,Lam)
    q = (s-u2);
    q_tilde=U'*q;
    r_tilde = q_tilde(1:end-2);    

    p = 2*t'*Phir*t/rho2/max(diag(Lam));

    q1 = q_tilde(end-1);
    q2 = q_tilde(end);           
    temp = 1+q2^2/q1^2;

    coeff =[1  -q1 0 0 -p/temp^2 ];    
    ro = roots(coeff);  % find roots for a polynomial with coefficient "coeff"

    for kk = length(ro):-1:1
       if abs(imag(ro(kk))-0)>=1e-4
           ro(kk) =[];
       end
    end    % remove complex roots

    for k =1:length(ro)
       if sign(ro(k)) == sign(q1)
           r_tilde(end+1) = ro(k);
       end
    end
    r_tilde(end+1) = r_tilde(end)*q2/q1;
    r=U*r_tilde;

