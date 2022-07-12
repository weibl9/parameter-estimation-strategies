function yf = nonlinearODE(sp,ic,tf)
% tf: time index
% ic: initial condition, at tf(1)
% sp: structural parameter

options = odeset('RelTol',3e-14,'AbsTol',3e-14);
f = @(t,s) sp(1)*s^2 + sp(2)*s + sp(3);
[~, yf] = ode45(f,tf,ic,options);
    
end