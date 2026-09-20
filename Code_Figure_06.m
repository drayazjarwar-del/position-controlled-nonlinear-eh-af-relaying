%% ARCHIVAL FILE: Code_Figure_06.m
% This archival copy preserves the numerical model, parameters, and algorithms
% of the supplied final script. Only dependency/output filenames were standardized
% for the reproducibility package.

%% FIGURE 6: Near-optimal power-splitting width vs normalized source transmit power
% Analytical result from Eq. (30) + direct Monte-Carlo validation
%
% The near-optimal width is
%
%   W_delta^(chi)
%     = |{rho_m in Q_rho : 10log10(Pout(rho_m)/Pout*) <= delta}| / N_rho
%
% where Pout* is the minimum outage over the SAME N_rho-point rho grid.
%
% Thus:
%   larger W_delta -> less sensitivity to the exact power-splitting ratio
%   smaller W_delta -> more accurate rho control is required
%
% IMPORTANT:
% 1) Run Code_Figure_02.m first.
% 2) This script loads Master_Numerical_Data.mat.
% 3) It therefore uses EXACTLY the same MA/PA realization and optimized
%    deterministic position gains as Figs. 2-5.
% 4) The analytical curves evaluate Eq. (30) over the complete
%    N_rho = 1001 grid at every normalized-source-transmit-power point.
% 5) Monte-Carlo markers use 10^6 common fading realizations and directly
%    evaluate Eqs. (7)-(10). To avoid an impractical 1001 x 10^6 matrix,
%    the empirical near-optimal boundaries are refined locally around the
%    analytical boundaries. The code first verifies that the analytical
%    near-optimal set is contiguous; otherwise it stops rather than making
%    an unsupported width approximation.
%
% FIGURE STYLE:
% - same one-column IEEE style as revised Fig. 2
% - main line width = 0.90
% - open simulation circles
% - light grid
% - 3.50 in figure width

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
    'snrDb', ...
    'Omega1_MA','Omega2_MA', ...
    'Omega1_PA','Omega2_PA', ...
    'M1','M2'};

for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k},'var')
        error('Required variable "%s" is missing from %s.', ...
              requiredVars{k},masterFile);
    end
end

%% ------------------------------------------------------------------------
% 2) TABLE-II / MANUSCRIPT PARAMETERS
% -------------------------------------------------------------------------
sigmaR2 = 1;
sigmaD2 = 1;

Pmax = sigmaD2*10^(13/10);      % Pmax/sigma_D^2 = 13 dB
bEH  = sigmaR2*10^(10/10);      % bEH/sigma_R^2 = 10 dB
aEH  = 0.5/sigmaR2;             % aEH*sigma_R^2 = 0.5
xi0  = 1/(1+exp(aEH*bEH));

R0 = 1;
gammaTh = 2^(2*R0)-1;

rhoMin = 0.01;
Nrho = 1001;
rhoGrid = linspace(rhoMin,1-rhoMin,Nrho);

% Near-optimality tolerances from Table II
deltaList = [0.1 0.5];          % dB

Nmc = 1e6;

% Default Log-alpha-mu tuple
alpha1 = 2; beta1 = 1; mu1 = 2; s1 = 0.5;
alpha2 = 2; beta2 = 1; mu2 = 2; s2 = 0.5;

% Adaptive Gauss-Legendre settings for the full outage-vs-rho curve
GLorders = [32 64 128 256 512];
relTol = 1e-9;
absTol = 1e-12;

p.alpha1=alpha1; p.beta1=beta1; p.mu1=mu1; p.s1=s1; p.M1=M1;
p.alpha2=alpha2; p.beta2=beta2; p.mu2=mu2; p.s2=s2; p.M2=M2;

p.sigmaR2=sigmaR2;
p.sigmaD2=sigmaD2;

p.Pmax=Pmax;
p.bEH=bEH;
p.aEH=aEH;
p.xi0=xi0;
p.gammaTh=gammaTh;

%% ------------------------------------------------------------------------
% 3) ANALYTICAL W_delta^(chi)
% -------------------------------------------------------------------------
nSNR = numel(snrDb);
nDelta = numel(deltaList);

Wana_MA = zeros(nDelta,nSNR);
Wana_PA = zeros(nDelta,nSNR);

rhoOpt_MA = zeros(1,nSNR);
rhoOpt_PA = zeros(1,nSNR);

Popt_MA = zeros(1,nSNR);
Popt_PA = zeros(1,nSNR);

% Store analytical boundary indices for Monte-Carlo refinement.
idxL_MA = zeros(nDelta,nSNR);
idxR_MA = zeros(nDelta,nSNR);
idxL_PA = zeros(nDelta,nSNR);
idxR_PA = zeros(nDelta,nSNR);

Nused_MA = zeros(1,nSNR);
Nused_PA = zeros(1,nSNR);

fprintf('\n===== FIG. 6: ANALYTICAL NEAR-OPTIMAL WIDTH =====\n');

for k = 1:nSNR

    thisSNR = snrDb(k);

    [PvecMA,Nused_MA(k)] = adaptiveOutageCurve( ...
        thisSNR,rhoGrid,Omega1_MA,Omega2_MA, ...
        p,GLorders,relTol,absTol);

    [PvecPA,Nused_PA(k)] = adaptiveOutageCurve( ...
        thisSNR,rhoGrid,Omega1_PA,Omega2_PA, ...
        p,GLorders,relTol,absTol);

    [Popt_MA(k),idxOptMA] = min(PvecMA);
    [Popt_PA(k),idxOptPA] = min(PvecPA);

    rhoOpt_MA(k) = rhoGrid(idxOptMA);
    rhoOpt_PA(k) = rhoGrid(idxOptPA);

    for d = 1:nDelta

        deltaDb = deltaList(d);

        % "Within delta dB of optimum" means:
        % 10log10(Pout/Popt) <= deltaDb.
        maskMA = 10*log10(PvecMA/Popt_MA(k)) <= deltaDb + 1e-12;
        maskPA = 10*log10(PvecPA/Popt_PA(k)) <= deltaDb + 1e-12;

        assertContiguous(maskMA, ...
            sprintf('MA, SNR=%g dB, delta=%g dB',thisSNR,deltaDb));

        assertContiguous(maskPA, ...
            sprintf('PA, SNR=%g dB, delta=%g dB',thisSNR,deltaDb));

        idMA = find(maskMA);
        idPA = find(maskPA);

        idxL_MA(d,k) = idMA(1);
        idxR_MA(d,k) = idMA(end);

        idxL_PA(d,k) = idPA(1);
        idxR_PA(d,k) = idPA(end);

        Wana_MA(d,k) = numel(idMA)/Nrho;
        Wana_PA(d,k) = numel(idPA)/Nrho;
    end

    fprintf(['SNR=%2d dB | rho*_MA=%.3f | ', ...
             'W_MA(0.1)=%.4f W_MA(0.5)=%.4f | ', ...
             'rho*_PA=%.3f | W_PA(0.1)=%.4f W_PA(0.5)=%.4f\n'], ...
             thisSNR,rhoOpt_MA(k), ...
             Wana_MA(1,k),Wana_MA(2,k), ...
             rhoOpt_PA(k), ...
             Wana_PA(1,k),Wana_PA(2,k));
end

%% ------------------------------------------------------------------------
% 4) MONTE-CARLO FADING REALIZATIONS
% -------------------------------------------------------------------------
rng(20260910,'twister');

r1 = @(u) ((1+beta1*u/mu1).^(1/beta1)-s1).^(1/alpha1) ...
          -(1-s1)^(1/alpha1);

r2 = @(u) ((1+beta2*u/mu2).^(1/beta2)-s2).^(1/alpha2) ...
          -(1-s2)^(1/alpha2);

U1mc = randg(mu1,Nmc,1);
U2mc = randg(mu2,Nmc,1);

Z1mc = r1(U1mc).^2/M1;
Z2mc = r2(U2mc).^2/M2;

%% ------------------------------------------------------------------------
% 5) MONTE-CARLO VALIDATION MARKERS
% -------------------------------------------------------------------------
% Use every 4 dB, matching the marker density adopted in revised Fig. 2.
simMask = mod(snrDb,4)==0;
simIdx = find(simMask);
snrSimDb = snrDb(simIdx);

Wsim_MA = nan(nDelta,nSNR);
Wsim_PA = nan(nDelta,nSNR);

fprintf('\n===== FIG. 6: MONTE-CARLO WIDTH VALIDATION =====\n');

for jj = 1:numel(simIdx)

    k = simIdx(jj);
    Ps = sigmaR2*10^(snrDb(k)/10);

    fprintf('SNR=%2d dB ...\n',snrDb(k));

    % Precompute effective fading gains for this architecture.
    X_MA = Omega1_MA*Z1mc;
    Y_MA = Omega2_MA*Z2mc;

    X_PA = Omega1_PA*Z1mc;
    Y_PA = Omega2_PA*Z2mc;

    for d = 1:nDelta

        deltaDb = deltaList(d);

        Wsim_MA(d,k) = empiricalNearOptimalWidthLocal( ...
            Ps,rhoGrid, ...
            idxL_MA(d,k),idxR_MA(d,k), ...
            rhoOpt_MA(k), ...
            deltaDb,X_MA,Y_MA,p);

        Wsim_PA(d,k) = empiricalNearOptimalWidthLocal( ...
            Ps,rhoGrid, ...
            idxL_PA(d,k),idxR_PA(d,k), ...
            rhoOpt_PA(k), ...
            deltaDb,X_PA,Y_PA,p);
    end

    fprintf(['  MA: W0.1=%.4f, W0.5=%.4f | ', ...
             'PA: W0.1=%.4f, W0.5=%.4f\n'], ...
             Wsim_MA(1,k),Wsim_MA(2,k), ...
             Wsim_PA(1,k),Wsim_PA(2,k));
end

%% ------------------------------------------------------------------------
% 6) CHECK WHERE THE WIDTH IS MINIMUM
% -------------------------------------------------------------------------
[~,kMinMA01] = min(Wana_MA(1,:));
[~,kMinMA05] = min(Wana_MA(2,:));
[~,kMinPA01] = min(Wana_PA(1,:));
[~,kMinPA05] = min(Wana_PA(2,:));

fprintf('\n===== MINIMUM NEAR-OPTIMAL WIDTH LOCATIONS =====\n');
fprintf('MA, delta=0.1 dB: %.1f dB\n',snrDb(kMinMA01));
fprintf('MA, delta=0.5 dB: %.1f dB\n',snrDb(kMinMA05));
fprintf('PA, delta=0.1 dB: %.1f dB\n',snrDb(kMinPA01));
fprintf('PA, delta=0.5 dB: %.1f dB\n',snrDb(kMinPA05));

%% ------------------------------------------------------------------------
% 7) REVISED IEEE-STYLE FIGURE 6
% -------------------------------------------------------------------------
cMA  = [1.00 0.00 0.00];
cPA  = [0.00 0.00 1.00];
cSIM = [0.00 0.00 0.00];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.50 2.55]);

ax = axes(fig);
hold(ax,'on');
box(ax,'on');

%% Analytical curves
% Same architecture colors as Fig. 2.
% Line style differentiates delta.

hMA01 = plot(ax, ...
    snrDb,Wana_MA(1,:),'-', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','MA, $\delta=0.1$ dB');

hMA05 = plot(ax, ...
    snrDb,Wana_MA(2,:),':', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','MA, $\delta=0.5$ dB');

hPA01 = plot(ax, ...
    snrDb,Wana_PA(1,:),'--', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','PA, $\delta=0.1$ dB');

hPA05 = plot(ax, ...
    snrDb,Wana_PA(2,:),'-.', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','PA, $\delta=0.5$ dB');

%% Simulation markers
% Open circles ensure the analytical lines remain visible through markers.
plot(ax, ...
    snrDb(simIdx),Wsim_MA(1,simIdx),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(simIdx),Wsim_MA(2,simIdx),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(simIdx),Wsim_PA(1,simIdx),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(simIdx),Wsim_PA(2,simIdx),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

%% Dummy simulation handle
hSim = plot(ax, ...
    NaN,NaN,'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'DisplayName','Simulation');

%% Axis labels
xlabel(ax, ...
    'Normalized source transmit power (dB)', ...
    'Interpreter','latex', ...
        'FontSize',10);

ylabel(ax, ...
    'Near-optimal width', ...
    'Interpreter','latex', ...
        'FontSize',10);

xlim(ax,[0 40]);
ylim(ax,[0 1]);

xticks(ax,0:5:40);
yticks(ax,0:0.2:1);

%% Clean IEEE-style axes/grid
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

%% Legend
lgd = legend(ax, ...
    [hMA01 hMA05 hPA01 hPA05 hSim], ...
    'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',7, ...
    'Box','on');

lgd.LineWidth = 0.55;
lgd.ItemTokenSize = [18 8];

%% Export
set(fig,'Renderer','painters');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_06.pdf'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_06.png'), ...
    'Resolution',600);

%% ------------------------------------------------------------------------
% 8) APPEND FIGURE-6 DATA TO MASTER FILE
% -------------------------------------------------------------------------
save(masterFile, ...
    'deltaList', ...
    'Wana_MA','Wana_PA', ...
    'Wsim_MA','Wsim_PA', ...
    'rhoOpt_MA','rhoOpt_PA', ...
    'Popt_MA','Popt_PA', ...
    'idxL_MA','idxR_MA','idxL_PA','idxR_PA', ...
    'Nused_MA','Nused_PA', ...
    '-append');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [PvecBest,Nused] = adaptiveOutageCurve( ...
    snrDb,rhoGrid,Omega1,Omega2,p,orders,relTol,absTol)
% Adaptive finite-order evaluation of Eq. (30) over the COMPLETE rho grid.

prev = [];
PvecBest = [];
Nused = orders(end);

for ii = 1:numel(orders)

    N = orders(ii);
    [xGL,wGL] = gaussLegendre(N);

    current = outageGLgrid( ...
        snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL);

    if ii > 1

        err = max(abs(current-prev));

        % Scale the convergence tolerance by the largest current outage.
        tol = absTol + relTol*max(abs(current));

        if err <= tol
            PvecBest = current;
            Nused = N;
            return;
        end
    end

    prev = current;
    PvecBest = current;
end

warning(['Full outage-vs-rho curve did not reach the requested ', ...
         'GL tolerance by N=%d at SNR %.1f dB. ', ...
         'Highest-order result is used.'], ...
         orders(end),snrDb);

end

%% ------------------------------------------------------------------------

function Pout = outageGLgrid( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL)
% Finite-order evaluation of Eq. (30) over all rho values simultaneously.

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

% Transform standard GL nodes from [-1,1] to [pth,1].
v = (1+pa)/2 + ((1-pa)/2).*xGL(:).';

v = min(v,1-eps);
v = max(v,realmin);

u = gammaincinv(v,p.mu1,'lower');

r = ((1+p.beta1*u/p.mu1).^(1/p.beta1)-p.s1) ...
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

Pout = min(max(Pout,realmin),1);

end

%% ------------------------------------------------------------------------

function W = empiricalNearOptimalWidthLocal( ...
    Ps,rhoGrid,idxLana,idxRana,rhoOptAna,deltaDb,X,Y,p)
% Monte-Carlo validation of W_delta using direct Eqs. (7)-(10).
%
% The analytical near-optimal set is first verified to be contiguous by the
% calling routine.  We then use its left/right boundaries as STARTING
% locations only, and move outward/inward using EMPIRICAL outage values
% until the empirical delta-dB crossings are found.
%
% This avoids evaluating 1001 million-sample outages at every SNR while
% still determining the empirical width from direct Monte Carlo samples.

Nrho = numel(rhoGrid);

% Small cache so repeated rho evaluations are not recomputed.
cache = nan(1,Nrho);

    function val = mcAt(idx)
        if isnan(cache(idx))
            cache(idx) = directMonteCarloOutageXY( ...
                Ps,rhoGrid(idx),X,Y,p);
        end
        val = cache(idx);
    end

% Empirical optimum:
% evaluate the analytical optimum and a small local neighborhood.
[~,idxOptAna] = min(abs(rhoGrid-rhoOptAna));

localOptIdx = max(1,idxOptAna-5):min(Nrho,idxOptAna+5);

localVals = zeros(size(localOptIdx));

for q = 1:numel(localOptIdx)
    localVals(q) = mcAt(localOptIdx(q));
end

[PmcBest,qMin] = min(localVals);
idxOptMC = localOptIdx(qMin);

threshold = PmcBest*10^(deltaDb/10);

% ---------------- Left empirical boundary ----------------
idx = min(idxLana,idxOptMC);

if mcAt(idx) <= threshold
    % Move left until just outside the near-optimal set.
    while idx > 1 && mcAt(idx-1) <= threshold
        idx = idx-1;
    end
    idxL = idx;
else
    % Move right until entering the near-optimal set.
    while idx < idxOptMC && mcAt(idx) > threshold
        idx = idx+1;
    end
    idxL = idx;
end

% ---------------- Right empirical boundary ----------------
idx = max(idxRana,idxOptMC);

if mcAt(idx) <= threshold
    % Move right until just outside the near-optimal set.
    while idx < Nrho && mcAt(idx+1) <= threshold
        idx = idx+1;
    end
    idxR = idx;
else
    % Move left until entering the near-optimal set.
    while idx > idxOptMC && mcAt(idx) > threshold
        idx = idx-1;
    end
    idxR = idx;
end

if idxR < idxL
    error('Empirical near-optimal interval became invalid.');
end

W = (idxR-idxL+1)/Nrho;

end

%% ------------------------------------------------------------------------

function Pout = directMonteCarloOutageXY(Ps,rho,X,Y,p)
% Direct Eqs. (7)-(10), with X and Y already formed.

Pin = rho*Ps.*X;
PH = EHmap(Pin,p);

gamma1 = (1-rho)*Ps.*X/p.sigmaR2;
gamma2 = PH.*Y/p.sigmaD2;

gammaE2E = ...
    (gamma1.*gamma2)./(gamma1+gamma2+1);

Pout = mean(gammaE2E < p.gammaTh);

end

%% ------------------------------------------------------------------------

function assertContiguous(mask,label)
% Stop if the analytical near-optimal set is not one contiguous interval.
% This protects the local Monte-Carlo boundary-refinement method.

idx = find(mask);

if isempty(idx)
    error('Near-optimal set is empty for %s.',label);
end

if any(diff(idx)~=1)
    error(['Near-optimal set is noncontiguous for %s. ', ...
           'A full Monte-Carlo rho-grid evaluation is required.'],label);
end

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
