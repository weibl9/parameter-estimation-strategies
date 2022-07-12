clc
clear
close all

%% data generation
n = 4;
h  = 1.5;                   % h in {0.5,1.0,1.5}
ts = h*(1:n)';

sp = [0.75 -0.25]; 
ic = 0.35;
ys = linearODE(sp,ic,ts);

%% parameter estimation
spic.true = [sp ic];     % true parameters
opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');

% ls-based estimates 
spic.ls = fit_ls(ts,ys,1); % 1 for linear ode 

% SNLS algorithm
t0 = 2*ts(1)-ts(2);
ic0 = 2*ys(1)-ys(2);

icl = -inf; icu = inf; 
lossfcn = @(p)fit_snls(p,ts,ys,t0,1); % 1 for linear ode  
ice = lsqnonlin(lossfcn,ic0,icl,icu,opts);
spic.sls = [spsnls(1:2) ice];
ysls_fit = linearODE(spic.sls(1:2),spic.sls(3),[t0;ts]);
spic.sls(3) = ysls_fit(2);   % ic is [t1, y1] not [t01 y0]

% NLS algorithm
po = [spic.ls(1:2) ys(1)];
pl = [-inf -inf -inf];
pu = [ inf  inf  inf];      % bounds
lossfcn = @(p)linearODE(p(1:2),p(3),ts)-ys;
spic.nls = lsqnonlin(lossfcn,po,pl,pu,opts);
ynls_fit = linearODE(spic.nls(1:2),spic.nls(3),ts);

%%
format long g
cell2mat(struct2cell(spic))

tf = (ts(1):h/10:8)';              % full time 
y.tr = linearODE(spic.true(1:2),  spic.true(3),tf);
y.ls = linearODE(spic.ls(1:2),    spic.ls(3),  tf);
y.nls  = linearODE(spic.nls(1:2), spic.nls(3), tf);
y.sls  = linearODE(spic.sls(1:2), spic.sls(3), tf);

figure
plot(ts, ys,'.-k','markersize',20); hold on
plot(tf, y.tr,'-r','linewidth',2);            % full 
plot(tf, y.ls,':','color','#7E2F8E','linewidth',2);
plot(tf, y.sls,'-.g','linewidth',2);     % time instants
plot(tf, y.nls,'--b','linewidth',2); hold off
xlim([0 inf])
xlabel('$t$','interpreter','latex'); 
ylabel('$y(t)$','interpreter','latex'); grid on 
legend('Observations','True trajectory', '2LS  trajectory',...
    'SLS  trajectory','NLS  trajectory','location','northwest')
xline(ts(end),'-k',{'Forecasting','Horizon'},'LineWidth',1,'HandleVisibility','off');

%% 



