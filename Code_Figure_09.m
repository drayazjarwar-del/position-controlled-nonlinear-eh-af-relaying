%% ARCHIVAL FILE: Code_Figure_09.m
% This archival copy preserves the numerical model, parameters, and algorithms
% of the supplied final script. Only dependency/output filenames were standardized
% for the reproducibility package.

%% FIGURE 9: Optimality verification and outage-evaluation reduction
% Real analytical verification of the proposed decomposed design against
% FULL JOINT EXHAUSTIVE SEARCH using the SAME rho grid.
%
% Figure purpose:
%   1) Verify that the proposed decomposition gives exactly the same
%      optimized outage as full joint exhaustive search.
%   2) Show the reduction in outage-objective evaluations:
%
%         Proposed:   N_rho
%         Exhaustive: K^2 N_rho
%
%      after the deterministic source/destination gain searches.
%
%      Hence
%
%         N_eval,exh / N_eval,prop = K^2.
%
% IMPORTANT:
% 1) Run Code_Figure_02.m first.
% 2) This script loads the SAME MA/PA spatial realizations saved by Fig. 2.
% 3) Normalized source transmit power is fixed at P_S/sigma_R^2 = 20 dB, matching the current
%    Fig. 9 setup in the manuscript.
% 4) N_rho = 1001, exactly as in Table II and the manuscript.
% 5) "Full exhaustive search" here is genuine: for each K it explicitly
%    enumerates all K^2 source/destination candidate-position pairs and
%    evaluates every rho in the 1001-point grid.
% 6) Candidate sets are NESTED subsets of the 41-point master grids.
%    K=1 is the normalized fixed reference position.  As K increases,
%    new spatially separated candidates are added without removing the
%    earlier ones. This makes the K-sweep reproducible and ensures that
%    increasing K represents adding position freedom.
%
% FIGURE STYLE:
% - Same restrained IEEE visual language as revised Fig. 2.
% - MA: red solid
% - PA: blue dashed
% - Fixed positions: black dotted
% - Full exhaustive search: black cross markers on BOTH MA and PA curves
% - Complexity ratio: dark gray dash-dot on the right y-axis
% - Analytical line width = 0.90
% - 3.50-inch single-column figure
%
% NOTE:
% The exhaustive computation is substantially heavier than Figs. 2-7
% because it explicitly evaluates K^2*N_rho designs. Progress is printed.

clear; clc; close all;

%% ------------------------------------------------------------------------
% 1) LOAD MASTER DATA FROM FIGURE 2
% -------------------------------------------------------------------------
% Resolve this script's own folder so shared data are found reliably.
thisScript = mfilename('fullpath');
if isempty(thisScript)
    thisDir = pwd;
else
    thisDir = fileparts(thisScript);
end
masterFile = fullfile(thisDir,'Master_Numerical_Data.mat');

% If the shared master data do not yet exist, generate them automatically
% from the finalized Figure-2 program in the same folder.
if ~isfile(masterFile)
    fig2File = fullfile(thisDir,'Code_Figure_02.m');

    if ~isfile(fig2File)
        error(['Master_Numerical_Data.mat was not found, and ', ...
               'Code_Figure_02.m is not in the same folder. ', ...
               'Keep all finalized Figure 2-9 programs together.']);
    end

    disp('Shared numerical data were not found.');
    disp('Running finalized Figure 2 first to generate them...');
    run(fig2File);
    close all;

    % Re-establish this script's folder after returning from Figure 2.
    thisScript = mfilename('fullpath');
    if isempty(thisScript)
        thisDir = pwd;
    else
        thisDir = fileparts(thisScript);
    end
    masterFile = fullfile(thisDir,'Master_Numerical_Data.mat');

    if ~isfile(masterFile)
        error(['Figure 2 completed, but Master_Numerical_Data.mat ', ...
               'was not created in the script folder.']);
    end
end

load(masterFile);

requiredVars = { ...
    'OmegaMA1','OmegaMA2', ...
    'OmegaPA1','OmegaPA2', ...
    'qMA','qPA', ...
    'M1','M2'};

for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k},'var')
        error('Required variable "%s" is missing from %s.', ...
              requiredVars{k},masterFile);
    end
end

%% ------------------------------------------------------------------------
% 2) PARAMETERS
% -------------------------------------------------------------------------
sigmaR2 = 1;
sigmaD2 = 1;

snrTestDb = 20;
Ps = sigmaR2*10^(snrTestDb/10);

Pmax = sigmaD2*10^(13/10);
bEH  = sigmaR2*10^(10/10);
aEH  = 0.5/sigmaR2;
xi0  = 1/(1+exp(aEH*bEH));

R0 = 1;
gammaTh = 2^(2*R0)-1;

rhoMin = 0.01;
Nrho = 1001;
rhoGrid = linspace(rhoMin,1-rhoMin,Nrho);

% Default Log-alpha-mu tuple
alpha1 = 2; beta1 = 1; mu1 = 2; s1 = 0.5;
alpha2 = 2; beta2 = 1; mu2 = 2; s2 = 0.5;

p.alpha1=alpha1; p.beta1=beta1; p.mu1=mu1; p.s1=s1; p.M1=M1;
p.alpha2=alpha2; p.beta2=beta2; p.mu2=mu2; p.s2=s2; p.M2=M2;
p.sigmaR2=sigmaR2; p.sigmaD2=sigmaD2;
p.Pmax=Pmax; p.bEH=bEH; p.aEH=aEH; p.xi0=xi0;
p.gammaTh=gammaTh;

% Candidate-position counts in the current Fig. 9.
Klist = [1 4 8 12 16 20];
Kmax = max(Klist);

% Fixed quadrature order for the exhaustive enumeration.
% A higher-order independent check is carried out at the final optima.
NqExh = 128;
[xGL,wGL] = gaussLegendre(NqExh);

% High-order final verification
NqCheck = 256;
[xGLcheck,wGLcheck] = gaussLegendre(NqCheck);

%% ------------------------------------------------------------------------
% 3) BUILD REPRODUCIBLE NESTED CANDIDATE SETS
% -------------------------------------------------------------------------
NfullMA = numel(OmegaMA1);
NfullPA = numel(OmegaPA1);

if NfullMA < Kmax || NfullPA < Kmax
    error('The master spatial grids must contain at least Kmax=%d points.',Kmax);
end

% Greedy farthest-point ordering in INDEX space.
% It starts from index 1, i.e., the normalized fixed reference position.
orderMA = nestedSpaceFillingOrder(NfullMA,Kmax);
orderPA = nestedSpaceFillingOrder(NfullPA,Kmax);

fprintf('\n===== NESTED CANDIDATE-SET ORDER =====\n');
fprintf('MA master-grid indices used up to K=20:\n');
disp(orderMA);
fprintf('PA master-grid indices used up to K=20:\n');
disp(orderPA);

%% ------------------------------------------------------------------------
% 4) FIXED-POSITION BENCHMARK
% -------------------------------------------------------------------------
Omega1_FIX = 1;
Omega2_FIX = 1;

[Pfixed,rhoFixed] = optimizeRhoFixedOrder( ...
    snrTestDb,Omega1_FIX,Omega2_FIX,rhoGrid,p,xGL,wGL);

fprintf('\n===== FIXED BENCHMARK =====\n');
fprintf('Pout_fixed = %.10e, rho*_fixed = %.4f\n',Pfixed,rhoFixed);

%% ------------------------------------------------------------------------
% 5) PROPOSED DECOMPOSED DESIGN AND FULL EXHAUSTIVE SEARCH
% -------------------------------------------------------------------------
nK = numel(Klist);

Pprop_MA = zeros(1,nK);
Pprop_PA = zeros(1,nK);

Pexh_MA = zeros(1,nK);
Pexh_PA = zeros(1,nK);

rhoProp_MA = zeros(1,nK);
rhoProp_PA = zeros(1,nK);

rhoExh_MA = zeros(1,nK);
rhoExh_PA = zeros(1,nK);

bestSrcExh_MA = zeros(1,nK);
bestDstExh_MA = zeros(1,nK);

bestSrcExh_PA = zeros(1,nK);
bestDstExh_PA = zeros(1,nK);

bestSrcProp_MA = zeros(1,nK);
bestDstProp_MA = zeros(1,nK);

bestSrcProp_PA = zeros(1,nK);
bestDstProp_PA = zeros(1,nK);

fprintf('\n===== FIG. 9: PROPOSED VS FULL EXHAUSTIVE SEARCH =====\n');

for kk = 1:nK

    K = Klist(kk);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('K = %d candidate positions per side\n',K);
    fprintf('------------------------------------------------------------\n');

    %% ---------------- MA candidate subset ----------------
    idxMA = orderMA(1:K);

    g1MA = OmegaMA1(idxMA);
    g2MA = OmegaMA2(idxMA);

    % Proposed decomposition:
    % maximize source-side and destination-side deterministic gains
    % independently, then perform one scalar rho search.
    [Omega1PropMA,i1local] = max(g1MA);
    [Omega2PropMA,i2local] = max(g2MA);

    bestSrcProp_MA(kk) = idxMA(i1local);
    bestDstProp_MA(kk) = idxMA(i2local);

    [Pprop_MA(kk),rhoProp_MA(kk)] = optimizeRhoFixedOrder( ...
        snrTestDb,Omega1PropMA,Omega2PropMA, ...
        rhoGrid,p,xGL,wGL);

    fprintf('MA proposed: Pout=%.10e, rho*=%.4f\n', ...
        Pprop_MA(kk),rhoProp_MA(kk));

    % Genuine full joint exhaustive search:
    [Pexh_MA(kk),rhoExh_MA(kk),i1exhMA,i2exhMA] = ...
        fullJointExhaustive( ...
        snrTestDb,rhoGrid,g1MA,g2MA,p,xGL,wGL, ...
        sprintf('MA, K=%d',K));

    bestSrcExh_MA(kk) = idxMA(i1exhMA);
    bestDstExh_MA(kk) = idxMA(i2exhMA);

    fprintf(['MA exhaustive: Pout=%.10e, rho*=%.4f, ', ...
             '|Delta|=%.3e\n'], ...
             Pexh_MA(kk),rhoExh_MA(kk), ...
             abs(Pprop_MA(kk)-Pexh_MA(kk)));

    %% ---------------- PA candidate subset ----------------
    idxPA = orderPA(1:K);

    g1PA = OmegaPA1(idxPA);
    g2PA = OmegaPA2(idxPA);

    [Omega1PropPA,i1local] = max(g1PA);
    [Omega2PropPA,i2local] = max(g2PA);

    bestSrcProp_PA(kk) = idxPA(i1local);
    bestDstProp_PA(kk) = idxPA(i2local);

    [Pprop_PA(kk),rhoProp_PA(kk)] = optimizeRhoFixedOrder( ...
        snrTestDb,Omega1PropPA,Omega2PropPA, ...
        rhoGrid,p,xGL,wGL);

    fprintf('PA proposed: Pout=%.10e, rho*=%.4f\n', ...
        Pprop_PA(kk),rhoProp_PA(kk));

    [Pexh_PA(kk),rhoExh_PA(kk),i1exhPA,i2exhPA] = ...
        fullJointExhaustive( ...
        snrTestDb,rhoGrid,g1PA,g2PA,p,xGL,wGL, ...
        sprintf('PA, K=%d',K));

    bestSrcExh_PA(kk) = idxPA(i1exhPA);
    bestDstExh_PA(kk) = idxPA(i2exhPA);

    fprintf(['PA exhaustive: Pout=%.10e, rho*=%.4f, ', ...
             '|Delta|=%.3e\n'], ...
             Pexh_PA(kk),rhoExh_PA(kk), ...
             abs(Pprop_PA(kk)-Pexh_PA(kk)));
end

%% ------------------------------------------------------------------------
% 6) HIGH-ORDER CHECK AT ALL FINAL OPTIMA
% -------------------------------------------------------------------------
% This verifies that using NqExh=128 for the exhaustive enumeration did
% not materially change the selected minimum compared with Nq=256.

Pprop_MA_check = zeros(1,nK);
Pprop_PA_check = zeros(1,nK);
Pexh_MA_check  = zeros(1,nK);
Pexh_PA_check  = zeros(1,nK);

for kk = 1:nK

    K = Klist(kk);

    idxMA = orderMA(1:K);
    idxPA = orderPA(1:K);

    % Proposed selected gains
    g1 = OmegaMA1(bestSrcProp_MA(kk));
    g2 = OmegaMA2(bestDstProp_MA(kk));
    [Pprop_MA_check(kk),~] = optimizeRhoFixedOrder( ...
        snrTestDb,g1,g2,rhoGrid,p,xGLcheck,wGLcheck);

    g1 = OmegaPA1(bestSrcProp_PA(kk));
    g2 = OmegaPA2(bestDstProp_PA(kk));
    [Pprop_PA_check(kk),~] = optimizeRhoFixedOrder( ...
        snrTestDb,g1,g2,rhoGrid,p,xGLcheck,wGLcheck);

    % Exhaustive-selected gains
    g1 = OmegaMA1(bestSrcExh_MA(kk));
    g2 = OmegaMA2(bestDstExh_MA(kk));
    [Pexh_MA_check(kk),~] = optimizeRhoFixedOrder( ...
        snrTestDb,g1,g2,rhoGrid,p,xGLcheck,wGLcheck);

    g1 = OmegaPA1(bestSrcExh_PA(kk));
    g2 = OmegaPA2(bestDstExh_PA(kk));
    [Pexh_PA_check(kk),~] = optimizeRhoFixedOrder( ...
        snrTestDb,g1,g2,rhoGrid,p,xGLcheck,wGLcheck);
end

fprintf('\n===== Nq=256 FINAL-OPTIMUM CHECK =====\n');

fprintf('Max MA proposed-vs-exhaustive difference = %.4e\n', ...
    max(abs(Pprop_MA_check-Pexh_MA_check)));

fprintf('Max PA proposed-vs-exhaustive difference = %.4e\n', ...
    max(abs(Pprop_PA_check-Pexh_PA_check)));

fprintf('Max Nq refinement change (MA proposed) = %.4e\n', ...
    max(abs(Pprop_MA_check-Pprop_MA)));

fprintf('Max Nq refinement change (PA proposed) = %.4e\n', ...
    max(abs(Pprop_PA_check-Pprop_PA)));

% Use the high-order checked proposed values for the plotted analytical
% curves.  Exhaustive markers use their own high-order checked values.
Pplot_MA = Pprop_MA_check;
Pplot_PA = Pprop_PA_check;

PplotExh_MA = Pexh_MA_check;
PplotExh_PA = Pexh_PA_check;

%% ------------------------------------------------------------------------
% 7) OUTAGE-EVALUATION COUNTS AND RATIO
% -------------------------------------------------------------------------
% Following the manuscript complexity discussion, these counts concern
% OUTAGE-OBJECTIVE evaluations AFTER deterministic gain searches.

NevalProp = Nrho*ones(size(Klist));
NevalExh = (Klist.^2)*Nrho;

evalRatio = NevalExh./NevalProp;      % exactly K^2

fprintf('\n===== OUTAGE-EVALUATION REDUCTION =====\n');

for kk = 1:nK
    fprintf(['K=%2d | proposed=%7d | exhaustive=%9d | ', ...
             'ratio=%4.0f x\n'], ...
             Klist(kk),NevalProp(kk),NevalExh(kk),evalRatio(kk));
end

%% ------------------------------------------------------------------------
% 8) VERIFY PROPOSITION NUMERICALLY
% -------------------------------------------------------------------------
tolVerify = 1e-9;

maxDiffMA = max(abs(Pplot_MA-PplotExh_MA));
maxDiffPA = max(abs(Pplot_PA-PplotExh_PA));

fprintf('\n===== OPTIMALITY VERIFICATION =====\n');
fprintf('max |Pprop - Pexh|, MA = %.4e\n',maxDiffMA);
fprintf('max |Pprop - Pexh|, PA = %.4e\n',maxDiffPA);

if maxDiffMA > tolVerify || maxDiffPA > tolVerify
    warning(['Proposed/exhaustive values differ by more than %.1e. ', ...
             'Inspect the selected candidate indices and quadrature order.'], ...
             tolVerify);
else
    fprintf('Proposed and full exhaustive search coincide numerically.\n');
end

%% ------------------------------------------------------------------------
% 9) REVISED IEEE-STYLE FIGURE 9
% -------------------------------------------------------------------------
cMA  = [1.00 0.00 0.00];
cPA  = [0.00 0.00 1.00];
cFIX = [0.00 0.00 0.00];
cEXH  = [0.00 0.00 0.00];
cCOMP = [0.24 0.24 0.24];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.50 2.62]);

ax = axes(fig);
hold(ax,'on');
box(ax,'on');

%% Left y-axis: optimized outage
yyaxis(ax,'left');

hMA = plot(ax, ...
    Klist,Pplot_MA,'-', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','Optimized MA');

hPA = plot(ax, ...
    Klist,Pplot_PA,'--', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','Optimized PA');

hFIX = plot(ax, ...
    Klist,Pfixed*ones(size(Klist)),':', ...
    'Color',cFIX, ...
    'LineWidth',0.90, ...
    'DisplayName','Fixed positions');

% Full exhaustive-search markers on BOTH architecture results.
% The same marker style is used because they represent the same method.
hEXH = plot(ax, ...
    Klist,PplotExh_MA,'x', ...
    'LineStyle','none', ...
    'MarkerSize',5.0, ...
    'MarkerEdgeColor',cEXH, ...
    'LineWidth',0.90, ...
    'DisplayName','Full exhaustive search');

plot(ax, ...
    Klist,PplotExh_PA,'x', ...
    'LineStyle','none', ...
    'MarkerSize',5.0, ...
    'MarkerEdgeColor',cEXH, ...
    'LineWidth',0.90, ...
    'HandleVisibility','off');

ylabel(ax, ...
    'Optimized outage probability', ...
    'Interpreter','latex', ...
        'FontSize',10);

% Let the actual data determine a compact range.
allP = [Pplot_MA Pplot_PA PplotExh_MA PplotExh_PA Pfixed];

pMin = min(allP);
pMax = max(allP);

padP = 0.08*(pMax-pMin);

if padP <= 0
    padP = 0.01;
end

ylim(ax,[max(0,pMin-padP), min(1,pMax+padP)]);

%% Right y-axis: outage-evaluation ratio
yyaxis(ax,'right');

hRatio = plot(ax, ...
    Klist,evalRatio,'-.', ...
    'Color',cCOMP, ...
    'LineWidth',0.90, ...
    'DisplayName', ...
    'Outage-evaluation ratio');

ylabel(ax,'Outage-evaluation ratio','FontSize',10);

ylim(ax,[0 1.05*max(evalRatio)]);
yticks(ax,0:100:ceil(max(evalRatio)/100)*100);

%% Common x-axis
xlabel(ax, ...
    'Candidate positions per side', ...
    'Interpreter','latex', ...
        'FontSize',10);

xlim(ax,[1 20]);
xticks(ax,Klist);

%% IEEE-style axes
set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',8, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',0.70, ...
    'TickDir','in', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'XMinorGrid','off', ...
    'YMinorGrid','off', ...
    'GridAlpha',0.14, ...
    'Layer','top');

% Make both y-axis colors neutral rather than MATLAB's default blue/red.
ax.YAxis(1).Color = [0.15 0.15 0.15];
ax.YAxis(2).Color = [0.15 0.15 0.15];

%% Return to left axis before legend
yyaxis(ax,'left');

%% Legend
lgd = legend(ax, ...
    [hMA hPA hFIX hEXH hRatio], ...
    'Location','northeast', ...
    'Interpreter','latex', ...
    'FontSize',7, ...
    'Box','on');

lgd.LineWidth = 0.55;
lgd.ItemTokenSize = [16 8];

%% Export
set(fig,'Renderer','painters');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_09.pdf'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_09.png'), ...
    'Resolution',600);

%% ------------------------------------------------------------------------
% 10) APPEND FIGURE-9 RESULTS TO MASTER FILE
% -------------------------------------------------------------------------
save(masterFile, ...
    'Klist','snrTestDb','Nrho', ...
    'orderMA','orderPA', ...
    'Pfixed','rhoFixed', ...
    'Pprop_MA','Pprop_PA', ...
    'Pexh_MA','Pexh_PA', ...
    'Pprop_MA_check','Pprop_PA_check', ...
    'Pexh_MA_check','Pexh_PA_check', ...
    'rhoProp_MA','rhoProp_PA', ...
    'rhoExh_MA','rhoExh_PA', ...
    'bestSrcProp_MA','bestDstProp_MA', ...
    'bestSrcProp_PA','bestDstProp_PA', ...
    'bestSrcExh_MA','bestDstExh_MA', ...
    'bestSrcExh_PA','bestDstExh_PA', ...
    'NevalProp','NevalExh','evalRatio', ...
    'NqExh','NqCheck', ...
    '-append');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function order = nestedSpaceFillingOrder(N,Kmax)
% Deterministic nested candidate ordering.
%
% Start from index 1 (the normalized fixed benchmark), then repeatedly add
% the grid point whose minimum distance to the already-selected set is
% largest. Ties are resolved by the smallest index.
%
% This yields nested, progressively space-filling subsets.

selected = false(1,N);

order = zeros(1,Kmax);
order(1) = 1;
selected(1) = true;

grid = 1:N;

for k = 2:Kmax

    chosen = order(1:k-1);

    dmin = inf(1,N);

    for j = 1:numel(chosen)
        dmin = min(dmin,abs(grid-chosen(j)));
    end

    dmin(selected) = -Inf;

    [~,idx] = max(dmin);

    order(k) = idx;
    selected(idx) = true;
end

end

%% ------------------------------------------------------------------------

function [Pbest,rhoBest] = optimizeRhoFixedOrder( ...
    snrDb,Omega1,Omega2,rhoGrid,p,xGL,wGL)
% One-dimensional rho optimization for fixed deterministic gains.

Pvec = outageGLgrid( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL);

[Pbest,idx] = min(Pvec);
rhoBest = rhoGrid(idx);

end

%% ------------------------------------------------------------------------

function [Pbest,rhoBest,iBest,jBest] = fullJointExhaustive( ...
    snrDb,rhoGrid,g1List,g2List,p,xGL,wGL,label)
% TRUE full joint exhaustive search:
%
%   min over source position i,
%            destination position j,
%            rho in the SAME N_rho grid.
%
% The routine explicitly enumerates all K^2 position pairs.  For speed,
% all destination gains are evaluated in one 3-D batch for each source
% position, but no position pair is omitted.

Ps = p.sigmaR2*10^(snrDb/10);

K1 = numel(g1List);
K2 = numel(g2List);
Nrho = numel(rhoGrid);
Nq = numel(xGL);

Pbest = Inf;
rhoBest = NaN;
iBest = NaN;
jBest = NaN;

rho = rhoGrid(:);

% 1 x 1 x K2 destination gains
g2 = reshape(g2List,1,1,K2);

% Quadrature weights: 1 x Nq x 1
w3 = reshape(wGL(:),1,Nq,1);

for i = 1:K1

    Omega1 = g1List(i);

    kappa = (1-rho)*Ps/p.sigmaR2;
    xth = p.gammaTh./kappa;

    uth = p.mu1*PsiFun( ...
        sqrt(p.M1*xth/Omega1), ...
        p.alpha1,p.beta1,p.s1);

    pth = gammainc(uth,p.mu1,'lower');

    % Nrho x Nq
    v = (1+pth)/2 + ((1-pth)/2).*xGL(:).';

    v = min(v,1-eps);
    v = max(v,realmin);

    u = gammaincinv(v,p.mu1,'lower');

    r = ...
        ((1+p.beta1*u/p.mu1).^(1/p.beta1)-p.s1) ...
        .^(1/p.alpha1) ...
        -(1-p.s1)^(1/p.alpha1);

    x = Omega1*r.^2/p.M1;

    Pin = rho.*Ps.*x;
    PH = EHmap(Pin,p);

    nu = PH/p.sigmaD2;

    den = nu.*(kappa.*x-p.gammaTh);

    T = inf(size(x));
    good = den > 0;

    num = p.gammaTh*(kappa.*x+1);
    T(good) = num(good)./den(good);

    % Nrho x Nq x K2
    arg = reshape(T,Nrho,Nq,1)./g2;

    FY = FZ( ...
        arg, ...
        p.alpha2,p.beta2,p.mu2,p.s2,p.M2);

    % Nrho x K2
    integralPart = squeeze(sum(FY.*w3,2));

    if K2 == 1
        integralPart = integralPart(:);
    end

    current = ...
        pth + ((1-pth)/2).*integralPart;

    current = min(max(current,0),1);

    % Minimum over rho for every destination candidate.
    [pairMin,rhoIdx] = min(current,[],1);

    [srcBest,jLocal] = min(pairMin);

    if srcBest < Pbest

        Pbest = srcBest;
        iBest = i;
        jBest = jLocal;
        rhoBest = rhoGrid(rhoIdx(jLocal));
    end

    fprintf('  %s: source candidate %d/%d complete\n', ...
        label,i,K1);
end

end

%% ------------------------------------------------------------------------

function Pout = outageGLgrid( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL)
% Finite-order evaluation of Eq. (30) over all rho values.

Ps = p.sigmaR2*10^(snrDb/10);

rho = rhoGrid(:);
kappa = (1-rho)*Ps/p.sigmaR2;

xth = p.gammaTh./kappa;

uth = p.mu1*PsiFun( ...
    sqrt(p.M1*xth/Omega1), ...
    p.alpha1,p.beta1,p.s1);

pth = gammainc(uth,p.mu1,'lower');

Pout = ones(size(rho));

active = find(pth < 1-1e-14);

if isempty(active)
    return;
end

pa = pth(active);
ka = kappa(active);

v = (1+pa)/2 + ((1-pa)/2).*xGL(:).';

v = min(v,1-eps);
v = max(v,realmin);

u = gammaincinv(v,p.mu1,'lower');

r = ...
    ((1+p.beta1*u/p.mu1).^(1/p.beta1)-p.s1) ...
    .^(1/p.alpha1) ...
    -(1-p.s1)^(1/p.alpha1);

x = Omega1*r.^2/p.M1;

Pin = rho(active).*Ps.*x;
PH = EHmap(Pin,p);

nu = PH/p.sigmaD2;

den = nu.*(ka.*x-p.gammaTh);

T = inf(size(x));
good = den > 0;

num = p.gammaTh*(ka.*x+1);
T(good) = num(good)./den(good);

FY = FZ( ...
    T/Omega2, ...
    p.alpha2,p.beta2,p.mu2,p.s2,p.M2);

integralPart = sum(FY.*wGL(:).',2);

Pout(active) = ...
    pa + ((1-pa)/2).*integralPart;

Pout = min(max(Pout,0),1);

end

%% ------------------------------------------------------------------------

function PH = EHmap(Pin,p)
% Eq. (7)

PH = ...
    (p.Pmax./(1+exp(-p.aEH*(Pin-p.bEH))) ...
    -p.Pmax*p.xi0)/(1-p.xi0);

PH = max(PH,0);

end

%% ------------------------------------------------------------------------

function F = FZ(z,alpha,beta,mu,s,M)
% Eq. (13)

z = max(z,0);

u = mu*PsiFun( ...
    sqrt(M*z), ...
    alpha,beta,s);

F = gammainc(u,mu,'lower');

end

%% ------------------------------------------------------------------------

function y = PsiFun(r,alpha,beta,s)
% Eq. (12)

c = (1-s)^(1/alpha);

y = ...
    (((r+c).^alpha+s).^beta-1)/beta;

end

%% ------------------------------------------------------------------------

function [x,w] = gaussLegendre(N)
% Standard N-point Gauss-Legendre nodes and weights on [-1,1].

k = (1:N-1)';

b = k./sqrt(4*k.^2-1);

J = diag(b,1)+diag(b,-1);

[V,D] = eig(J);

[x,idx] = sort(diag(D));

V = V(:,idx);

w = 2*(V(1,:)').^2;

end
