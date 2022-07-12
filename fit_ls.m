function spic = fit_ls(t,y,type)
% t:  time index, column vector
% y:  Cusum vector,  column vector
% type: 1 for linear and 2 for nonlinear case

dt = diff(t);
dy = diff(y);
x  = dy./dt;

n = length(dy);

switch type
    case 1  % f(y) = ay+b
        Theta = ones(n,2);
        Theta(:,1) = (y(1:end-1)+y(2:end))/2;        
        sp = Theta\x; 
        
        opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
        ic0 = y(1);
        icl = -inf; icu = inf; 
        lossfcn2 = @(ic)linearODE(sp,ic,t)-y;
        ic = lsqnonlin(lossfcn2,ic0,icl,icu,opts);        
    case 2  % f(y) = ay^2+by+c
        Theta = ones(n,3);
        Theta(:,1) = (y(1:end-1).^2 + y(2:end).^2)/2; 
        Theta(:,2) = (y(1:end-1)    + y(2:end)   )/2;        
        sp = Theta\x;
        
        opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
        ic0 = y(1);
        icl = -inf; icu = inf; 
        lossfcn2 = @(ic)nonlinearODE(sp,ic,t)-y;
        ic = lsqnonlin(lossfcn2,ic0,icl,icu,opts); 
end

spic = [sp' ic];

end