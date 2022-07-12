clc
clear
close all

%%
nrep = 500;
rng('default'); 

tbl = [];

for nvr = [.15 .10 .05]
spicall = {nan(nrep,4), nan(nrep,4), nan(nrep,4)};  % 2ls, sls, nls

for irep = 1:nrep

    %% data generation
    n = 20;
    h  = 0.2;
    ts = h*(1:n)';              % full time 

    sp = [-0.8 1.2 0.1]; 
    ic = 0.1;
    yt = nonlinearODE(sp,ic,ts);

%     nvr = 0.1;
    rng(irep);
    ys = yt + nvr*std(yt)*randn(size(yt));

    %% parameter estimation
    spic.true = [sp ic];     % true parameters 
    opts = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');

    % ls-based estimates 
    spic.ls = fit_ls(ts,ys,2);  % 2 for nonlinear ode

    % SNLS algorithm
    t0 = 2*ts(1)-ts(2);
    ic0 = 2*ys(1)-ys(2);

    icl = -inf; icu = inf; 
    lossfcn = @(p)fit_snls(p,ts,ys,t0,2); % 2 for nonlinear ode
    ice = lsqnonlin(lossfcn,ic0,icl,icu,opts);
    spic.sls = [spsnls(1:3) ice];
    ysls_fit = nonlinearODE(spic.sls(1:3),spic.sls(4),[t0;ts]);
    spic.sls(4) = ysls_fit(2);   % ic is [t1, y1] not [t01 y0]

    % NLS algorithm
    po = [spic.ls(1:3) ys(1)];
    pl = [-inf -inf -inf -inf];
    pu = [ inf  inf  inf  inf];      % bounds
    lossfcn = @(p)nonlinearODE(p(1:3),p(4),ts)-ys;
    spic.nls = lsqnonlin(lossfcn,po,pl,pu,opts);
    ynls_fit = nonlinearODE(spic.nls(1:3),spic.nls(4),ts);
    
    %% save parameter estimates 
    spicall{1}(irep,:) = spic.ls; 
    spicall{2}(irep,:) = spic.sls; 
    spicall{3}(irep,:) = spic.nls; 

end

tbl = [tbl 
       [n  nvr mean(spicall{1,1}) mean(spicall{1,2}) mean(spicall{1,3}); 
       nan nan  std(spicall{1,1})  std(spicall{1,2})  std(spicall{1,3}) ] ];

end

tbl = tbl(:, [1 2 [3 7 11] [3 7 11]+1 [3 7 11]+2 [3 7 11]+3])

%%
