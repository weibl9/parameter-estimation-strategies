function yf = linearODE(sp,ic,tf)
% tf: time index
% ic: initial condition, at tf(1)
% sp: structural parameter

% syms y(t)
% eqn = diff(y,t) == sp(1)*y+sp(2);
% cond = y(tf(1)) == ic;
% yfcn = dsolve(eqn,cond);
% yf = eval(subs(yfcn,tf));

yf = -(sp(2) - sp(1)*exp(sp(1)*tf)*exp(-sp(1)*tf(1))*(ic + sp(2)/sp(1)))/sp(1);
    
end