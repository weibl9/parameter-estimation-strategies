function loss = fit_snls(ic,t,y,t0,type)
% ic: initial condition 
% t:  time index, column vector
% y:  Cusum vector,  column vector
% type: 1 for linear and 2 for nonlinear case

t_aug = [t0; t];
y_aug = [ic; y];

dt = diff(t_aug);
dy = diff(y_aug);

n  = length(dy);
x  = dy./dt;

switch type 
    case 1  % f(y) = ay+b
        Theta = ones(n,2);
        Theta(:,1) = (y_aug(1:end-1)+y_aug(2:end))/2;
        
        beta = Theta\x;
        loss = x - Theta*beta;
    case 2  % f(y) = ay^2+by+c
        Theta = ones(n,3);
        Theta(:,1) = (y_aug(1:end-1).^2 + y_aug(2:end).^2)/2;
        Theta(:,2) = (y_aug(1:end-1)    + y_aug(2:end)   )/2;
        
        beta = Theta\x;
        loss = x - Theta*beta;
end

assignin('base','spsnls',beta');

end




