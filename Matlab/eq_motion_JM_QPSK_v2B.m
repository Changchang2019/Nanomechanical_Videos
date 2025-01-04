%solve m ddz + Gamma dz+k z + alpha z^3 + eta z^2 dz= Fd(t)
function [t,x] = eq_motion_JM_QPSK_v2B(x0,time_span,omega0,Q,meff,alpha,eta,Famp,Fwave,dt,ode)
    % options = odeset('RelTol',1e-13,'AbsTol',[1e-13 1e-13]);
    % options = odeset('RelTol',1e-5);
    % [t,q] = ode45(@eqm,time_span,x0);
    switch ode
        case 23
        [t,x] = ode23(@eqm,time_span,x0);
        otherwise
        [t,x] = ode45(@eqm,time_span,x0);
    end

        
    function xdot=eqm(t,x)
    xdot_1 = x(2);
    xdot_2 = -omega0/Q*x(2)-omega0^2*x(1)+force(t)/meff- alpha/meff*x(1)^3-eta/meff*x(2)*x(1)^2;
    xdot = [xdot_1 ; xdot_2];
    end

    function f = force(t)
    %f = forceAmp*sin(omega_d*t);
    if 1+floor(t/dt)<= size(Fwave,1)
    f = Famp * Fwave(1+floor(t/dt));
    else
        f = 0;
    end
    end
end

