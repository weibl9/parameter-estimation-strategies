clc
clear
close all

%%
case_number = 2;

switch case_number
    case 1      %% LINEAR 
        %% data
        xall = [9.4, 12.5, 14.0, 15.9, 19.3, 24.1, 25.8, 28.7, 39.6, 42.2, 58.3, 77.5, 89.6, 98, 106.4]';
        ns = 6;
        xs = xall(1:ns); 
        xf = xall(ns+1:end);

        h = 1; 
        ts = h*[1:length(xs)]';
        tf = [ts; ts(end)+h*(1:length(xf))'];

        ys = cumsum(xs)*h;
        yall = cumsum(xall)*h;

        %% parameter estimation
        opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');

        spic.ls = fit_ls(ts,ys,1); % 1 for linear ode  

        % NLS algorithm
        po = [spic.ls(1:2) ys(1)];
        pl = [-inf -inf -inf];
        pu = [ inf  inf  inf];      % bounds
        lossfcn = @(p)linearODE(p(1:2),p(3),ts)-ys;
        spic.nls = lsqnonlin(lossfcn,po,pl,pu,opts);
        ynls_fit = linearODE(spic.nls(1:2),spic.nls(3),ts);

        % SNLS algorithm
        t0 = 2*ts(1)-ts(2);
        ic0 = 2*ys(1)-ys(2);

        icl = -inf; icu = inf; 
        lossfcn = @(p)fit_snls(p,ts,ys,t0,1); % 1 for linear ode  
        ice = lsqnonlin(lossfcn,ic0,icl,icu,opts);
        spic.sls = [spsnls(1:2) ice];
        ysls_fit = linearODE(spic.sls(1:2),spic.sls(3),[t0;ts]);
        spic.sls(3) = ysls_fit(2);   % ic is [t1, y1] not [t01 y0]

        %%
        format long g
        cell2mat(struct2cell(spic))

        y.fir = linearODE(spic.ls(1:2),  xall(1),  tf);
        y.ls  = linearODE(spic.ls(1:2),  spic.ls(3),  tf);
        y.nls = linearODE(spic.nls(1:2), spic.nls(3), tf);
        y.sls = linearODE(spic.sls(1:2), spic.sls(3), tf);    

        x.fir = [y.fir(1); diff(y.fir)]/h;
        x.ls  = [y.ls(1);  diff(y.ls) ]/h;
        x.nls = [y.nls(1); diff(y.nls)]/h;
        x.sls = [y.sls(1); diff(y.sls)]/h;

    case 2    % NONLINEAR
        xall = 0.01*[340 346 412 409 458 405 387 nan(1,8)]'; ns = 6;
        xs = xall(1:ns); 
        xf = xall(ns+1:end);

        h = 1; 
        ts = h*[1:length(xs)]';
        tf = [ts; ts(end)+h*(1:length(xf))'];

        ys = cumsum(xs)*h;
        yall = cumsum(xall)*h;

        %% parameter estimation
        opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');

        spic.ls = fit_ls(ts,ys,2);  % 2 for nonlinear ode

        % NLS algorithm
        po = [spic.ls(1:3) ys(1)];
        pl = [-inf -inf -inf -inf];
        pu = [ inf  inf  inf  inf];      % bounds
        lossfcn = @(p)nonlinearODE(p(1:3),p(4),ts)-ys;
        spic.nls = lsqnonlin(lossfcn,po,pl,pu,opts);
        ynls_fit = nonlinearODE(spic.nls(1:3),spic.nls(4),ts);

        % SNLS algorithm
        t0 = 2*ts(1)-ts(2);
        ic0 = 2*ys(1)-ys(2);

        icl = -inf; icu = inf; 
        lossfcn = @(p)fit_snls(p,ts,ys,t0,2); % 2 for nonlinear ode
        ice = lsqnonlin(lossfcn,ic0,icl,icu,opts);
        spic.sls = [spsnls(1:3) ice];
        ysls_fit = nonlinearODE(spic.sls(1:3),spic.sls(4),[t0;ts]);
        spic.sls(4) = ysls_fit(2);   % ic is [t1, y1] not [t01 y0]

        %%
        format long g
        spic

        y.fir = nonlinearODE(spic.ls(1:3),  xall(1),  tf);
        y.ls  = nonlinearODE(spic.ls(1:3),  spic.ls(4),  tf);
        y.nls = nonlinearODE(spic.nls(1:3), spic.nls(4), tf);
        y.sls = nonlinearODE(spic.sls(1:3), spic.sls(4), tf);   

        x.fir = [y.fir(1); diff(y.fir)]/h;
        x.ls  = [y.ls(1);  diff(y.ls) ]/h;
        x.nls = [y.nls(1); diff(y.nls)]/h;
        x.sls = [y.sls(1); diff(y.sls)]/h;
end

%% 
f = figure;

format short 
yape = [abs(y.fir-yall)./yall, abs(y.ls-yall)./yall, ...
        abs(y.sls-yall)./yall, abs(y.nls-yall)./yall ]*100;

xape = [abs(x.fir-xall)./xall, abs(x.ls-xall)./xall, ...
        abs(x.sls-xall)./xall, abs(x.nls-xall)./xall ]*100;

% [mean(xape(1:ns,:)) mean(xape(ns+1:end,:))]
    
subplot(121)
plot(tf, yall,'.k','markersize',25); hold on
plot(tf, y.ls,'-r','linewidth',2);            % full 
plot(tf, y.sls,'-.g','linewidth',2);     % time instants
plot(tf, y.nls,'--b','linewidth',2); hold off
xlabel('$t$','interpreter','latex'); 
ylabel('$y(t)$','interpreter','latex'); grid on 
legend('Observations','2LS model','SLS model','NLS model','location','best')
xline(ts(end)+h/2,'--m',{'Forecasts'},'LineWidth',1,'HandleVisibility','off');
xlim([tf(1)-h/2, tf(end)+h/2])
set(gca,'fontsize',12)

subplot(122)
plot(tf, xall,'.-k','markersize',25); hold on
plot(tf, x.ls,'-r','linewidth',2);            % full 
plot(tf, x.sls,'-.g','linewidth',2);     % time instants
plot(tf, x.nls,'--b','linewidth',2); hold off
xlabel('$t$','interpreter','latex'); 
ylabel('$x(t)$','interpreter','latex'); grid on 
legend('Observations','2LS model','SLS model','NLS model','location','best')
xline(ts(end)+h/2,'--m',{'Forecasts'},'LineWidth',1,'HandleVisibility','off');
xlim([tf(1)-h/2, tf(end)+h/2])

set(gca,'fontsize',12)
set(gcf,'position',[100 200 1300 450])






