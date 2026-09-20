%% FIGURE 7: Hop-specific outage sensitivity to Log-alpha-mu clustering
% Analytical result from Eq. (30) + direct Monte-Carlo validation
%
% Fig. 7 studies the effect of reducing the clustering parameter mu_i
% from its default value 2 to 1 on ONE hop at a time, while keeping the
% other hop at mu = 2.
%
% Define
%
%   S_i^(chi) = 10 log10( Pout^(chi)(mu_i=1) /
%                          Pout^(chi)(mu_1=mu_2=2) )   [dB].
%
% Hence:
%   S_i > 0  -> reducing mu_i from 2 to 1 increases outage;
%   S_i < 0  -> reducing mu_i gives a small outage reduction at that SNR,
%               possible at low SNR because unit-mean fading CDFs cross.
%
% IMPORTANT:
% 1) Run JSAC_Figure02_FINAL.m first.
% 2) This script loads JSAC_Master_Numerical_Data.mat.
% 3) It therefore uses EXACTLY the same optimized MA/PA deterministic
%    position gains and baseline outage results as Fig. 2.
% 4) For EACH perturbed fading case, rho is reoptimized over the SAME
%    N_rho = 1001 grid.
% 5) When mu changes, the normalization constant M_i is recomputed so that
%    the residual fading power remains unit mean. This is essential.
%
% FIGURE STYLE:
% - Same one-column IEEE style as revised Fig. 2.
% - Main analytical lines: 0.90 pt.
% - Open simulation circles: the analytical line remains visible inside.
% - MA and PA use the same restrained colors as Fig. 2.
% - Line style differentiates first-hop and second-hop sensitivity.

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
masterFile = fullfile(thisDir,'JSAC_Master_Numerical_Data.mat');

% If the shared master data do not yet exist, generate them automatically
% from the finalized Figure-2 program in the same folder.
if ~isfile(masterFile)
    fig2File = fullfile(thisDir,'JSAC_Figure02_FINAL.m');

    if ~isfile(fig2File)
        error(['JSAC_Master_Numerical_Data.mat was not found, and ', ...
               'JSAC_Figure02_FINAL.m is not in the same folder. ', ...
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
    masterFile = fullfile(thisDir,'JSAC_Master_Numerical_Data.mat');

    if ~isfile(masterFile)
        error(['Figure 2 completed, but JSAC_Master_Numerical_Data.mat ', ...
               'was not created in the script folder.']);
    end
end

load(masterFile);

requiredVars = { ...
    'snrDb', ...
    'Pana_MA','Pana_PA', ...
    'Psim_MA','Psim_PA', ...
    'Omega1_MA','Omega2_MA', ...
    'Omega1_PA','Omega2_PA', ...
    'M1','M2'};

for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k},'var')
        error('Required variable "%s" is missing from %s.', ...
              requiredVars{k},masterFile);
    end
end

% Preserve the baseline normalization constants from Fig. 2.
M1_base = M1;
M2_base = M2;

%% ------------------------------------------------------------------------
% 2) TABLE-II / BASELINE PARAMETERS
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

Nmc = 1e6;

% Default Log-alpha-mu parameters
alpha1 = 2; beta1 = 1; mu1_base = 2; s1 = 0.5;
alpha2 = 2; beta2 = 1; mu2_base = 2; s2 = 0.5;

% Perturbed clustering parameter
muReduced = 1;

% Gauss-Legendre convergence settings
relTol = 1e-9;
absTol = 1e-12;
GLorders = [32 64 128 256 512];

%% ------------------------------------------------------------------------
% 3) RECOMPUTE UNIT-MEAN NORMALIZATION FOR mu = 1
% -------------------------------------------------------------------------
% Baseline M1_base and M2_base correspond to mu1=mu2=2 and were produced
% by the Fig. 2 script.  When mu_i changes, M_i must also change.

M1_mu1 = normalizationM(alpha1,beta1,muReduced,s1);
M2_mu1 = normalizationM(alpha2,beta2,muReduced,s2);

fprintf('\n===== UNIT-MEAN NORMALIZATION CONSTANTS =====\n');
fprintf('Hop 1: M(mu=2) = %.10f, M(mu=1) = %.10f\n', ...
    M1_base,M1_mu1);
fprintf('Hop 2: M(mu=2) = %.10f, M(mu=1) = %.10f\n', ...
    M2_base,M2_mu1);

%% ------------------------------------------------------------------------
% 4) PARAMETER STRUCTURES FOR THE TWO PERTURBED CASES
% -------------------------------------------------------------------------
% Case A: first hop changed, mu1 = 1; second hop remains mu2 = 2
pFirst.alpha1 = alpha1;
pFirst.beta1  = beta1;
pFirst.mu1    = muReduced;
pFirst.s1     = s1;
pFirst.M1     = M1_mu1;

pFirst.alpha2 = alpha2;
pFirst.beta2  = beta2;
pFirst.mu2    = mu2_base;
pFirst.s2     = s2;
pFirst.M2     = M2_base;

pFirst.sigmaR2 = sigmaR2;
pFirst.sigmaD2 = sigmaD2;
pFirst.Pmax = Pmax;
pFirst.bEH = bEH;
pFirst.aEH = aEH;
pFirst.xi0 = xi0;
pFirst.gammaTh = gammaTh;

% Case B: second hop changed, mu2 = 1; first hop remains mu1 = 2
pSecond.alpha1 = alpha1;
pSecond.beta1  = beta1;
pSecond.mu1    = mu1_base;
pSecond.s1     = s1;
pSecond.M1     = M1_base;

pSecond.alpha2 = alpha2;
pSecond.beta2  = beta2;
pSecond.mu2    = muReduced;
pSecond.s2     = s2;
pSecond.M2     = M2_mu1;

pSecond.sigmaR2 = sigmaR2;
pSecond.sigmaD2 = sigmaD2;
pSecond.Pmax = Pmax;
pSecond.bEH = bEH;
pSecond.aEH = aEH;
pSecond.xi0 = xi0;
pSecond.gammaTh = gammaTh;

%% ------------------------------------------------------------------------
% 5) ANALYTICAL OUTAGE WITH mu_i REDUCED TO 1
% -------------------------------------------------------------------------
nSNR = numel(snrDb);

% Perturbed analytical outage
Pana_MA_first  = zeros(1,nSNR);
Pana_MA_second = zeros(1,nSNR);
Pana_PA_first  = zeros(1,nSNR);
Pana_PA_second = zeros(1,nSNR);

% Reoptimized rho in each perturbed case
rho_MA_first  = zeros(1,nSNR);
rho_MA_second = zeros(1,nSNR);
rho_PA_first  = zeros(1,nSNR);
rho_PA_second = zeros(1,nSNR);

fprintf('\n===== FIG. 7: ANALYTICAL PERTURBED OUTAGE =====\n');

for k = 1:nSNR

    thisSNR = snrDb(k);

    [Pana_MA_first(k),rho_MA_first(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_MA,Omega2_MA, ...
        rhoGrid,pFirst,GLorders,relTol,absTol);

    [Pana_MA_second(k),rho_MA_second(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_MA,Omega2_MA, ...
        rhoGrid,pSecond,GLorders,relTol,absTol);

    [Pana_PA_first(k),rho_PA_first(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_PA,Omega2_PA, ...
        rhoGrid,pFirst,GLorders,relTol,absTol);

    [Pana_PA_second(k),rho_PA_second(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_PA,Omega2_PA, ...
        rhoGrid,pSecond,GLorders,relTol,absTol);

    fprintf(['SNR=%2d dB | ', ...
             'MA first %.4e, second %.4e | ', ...
             'PA first %.4e, second %.4e\n'], ...
             thisSNR, ...
             Pana_MA_first(k),Pana_MA_second(k), ...
             Pana_PA_first(k),Pana_PA_second(k));
end

%% ------------------------------------------------------------------------
% 6) ANALYTICAL HOP-SPECIFIC OUTAGE CHANGE
% -------------------------------------------------------------------------
% Baseline analytical Pana_MA / Pana_PA are loaded from Fig. 2 and use
% mu1=mu2=2 with rho optimized for that baseline case.

Sana_MA_first = ...
    10*log10(Pana_MA_first./Pana_MA);

Sana_MA_second = ...
    10*log10(Pana_MA_second./Pana_MA);

Sana_PA_first = ...
    10*log10(Pana_PA_first./Pana_PA);

Sana_PA_second = ...
    10*log10(Pana_PA_second./Pana_PA);

%% ------------------------------------------------------------------------
% 7) MONTE-CARLO FADING REALIZATIONS
% -------------------------------------------------------------------------
% Reproduce the baseline Fig. 2 residual-fading samples exactly.
rng(20260910,'twister');

r1_base = @(u) ...
    ((1+beta1*u/mu1_base).^(1/beta1)-s1).^(1/alpha1) ...
    -(1-s1)^(1/alpha1);

r2_base = @(u) ...
    ((1+beta2*u/mu2_base).^(1/beta2)-s2).^(1/alpha2) ...
    -(1-s2)^(1/alpha2);

U1_base = randg(mu1_base,Nmc,1);
U2_base = randg(mu2_base,Nmc,1);

Z1_base = r1_base(U1_base).^2/M1_base;
Z2_base = r2_base(U2_base).^2/M2_base;

% Perturbed first-hop samples (mu1 = 1).
rng(20260911,'twister');

r1_red = @(u) ...
    ((1+beta1*u/muReduced).^(1/beta1)-s1).^(1/alpha1) ...
    -(1-s1)^(1/alpha1);

U1_red = randg(muReduced,Nmc,1);
Z1_red = r1_red(U1_red).^2/M1_mu1;

% Perturbed second-hop samples (mu2 = 1).
rng(20260912,'twister');

r2_red = @(u) ...
    ((1+beta2*u/muReduced).^(1/beta2)-s2).^(1/alpha2) ...
    -(1-s2)^(1/alpha2);

U2_red = randg(muReduced,Nmc,1);
Z2_red = r2_red(U2_red).^2/M2_mu1;

%% ------------------------------------------------------------------------
% 8) MONTE-CARLO OUTAGE FOR EACH PERTURBED CASE
% -------------------------------------------------------------------------
Psim_MA_first  = zeros(1,nSNR);
Psim_MA_second = zeros(1,nSNR);
Psim_PA_first  = zeros(1,nSNR);
Psim_PA_second = zeros(1,nSNR);

% Also recompute the baseline MC values using exactly the same samples
% and the baseline rho values saved by Fig. 2. This makes the dB-ratio
% definition internally consistent inside this script.
PsimBase_MA = zeros(1,nSNR);
PsimBase_PA = zeros(1,nSNR);

% Baseline rho values must be present from Fig. 2.
if ~exist('rho_MA','var') || ~exist('rho_PA','var')
    error(['rho_MA/rho_PA are missing from the master file. ', ...
           'Run JSAC_Figure02_FINAL.m first.']);
end

fprintf('\n===== FIG. 7: MONTE-CARLO VALIDATION =====\n');

for k = 1:nSNR

    Ps = sigmaR2*10^(snrDb(k)/10);

    % Baseline mu1=mu2=2
    PsimBase_MA(k) = directMonteCarloOutage( ...
        Ps,rho_MA(k),Omega1_MA,Omega2_MA, ...
        Z1_base,Z2_base,pFirst);  % common system params are identical

    PsimBase_PA(k) = directMonteCarloOutage( ...
        Ps,rho_PA(k),Omega1_PA,Omega2_PA, ...
        Z1_base,Z2_base,pFirst);

    % First-hop perturbation: mu1=1, mu2=2
    Psim_MA_first(k) = directMonteCarloOutage( ...
        Ps,rho_MA_first(k),Omega1_MA,Omega2_MA, ...
        Z1_red,Z2_base,pFirst);

    Psim_PA_first(k) = directMonteCarloOutage( ...
        Ps,rho_PA_first(k),Omega1_PA,Omega2_PA, ...
        Z1_red,Z2_base,pFirst);

    % Second-hop perturbation: mu1=2, mu2=1
    Psim_MA_second(k) = directMonteCarloOutage( ...
        Ps,rho_MA_second(k),Omega1_MA,Omega2_MA, ...
        Z1_base,Z2_red,pSecond);

    Psim_PA_second(k) = directMonteCarloOutage( ...
        Ps,rho_PA_second(k),Omega1_PA,Omega2_PA, ...
        Z1_base,Z2_red,pSecond);

    fprintf(['SNR=%2d dB | ', ...
             'MA first %.4e, second %.4e | ', ...
             'PA first %.4e, second %.4e\n'], ...
             snrDb(k), ...
             Psim_MA_first(k),Psim_MA_second(k), ...
             Psim_PA_first(k),Psim_PA_second(k));
end

%% ------------------------------------------------------------------------
% 9) MONTE-CARLO HOP-SPECIFIC OUTAGE CHANGE
% -------------------------------------------------------------------------
Ssim_MA_first = ...
    10*log10(Psim_MA_first./PsimBase_MA);

Ssim_MA_second = ...
    10*log10(Psim_MA_second./PsimBase_MA);

Ssim_PA_first = ...
    10*log10(Psim_PA_first./PsimBase_PA);

Ssim_PA_second = ...
    10*log10(Psim_PA_second./PsimBase_PA);

%% ------------------------------------------------------------------------
% 10) HIGH-SOURCE-POWER LIMITS FROM EQ. (34)
% -------------------------------------------------------------------------
% First-hop fading disappears from Eq. (34), hence:
Slim_MA_first = 0;
Slim_PA_first = 0;

% For the second hop:
%   S2_inf = 10 log10[Pfloor(mu2=1)/Pfloor(mu2=2)].

Pfloor_MA_mu2base = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_MA), ...
    alpha2,beta2,mu2_base,s2,M2_base);

Pfloor_MA_mu2red = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_MA), ...
    alpha2,beta2,muReduced,s2,M2_mu1);

Pfloor_PA_mu2base = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_PA), ...
    alpha2,beta2,mu2_base,s2,M2_base);

Pfloor_PA_mu2red = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_PA), ...
    alpha2,beta2,muReduced,s2,M2_mu1);

Slim_MA_second = ...
    10*log10(Pfloor_MA_mu2red/Pfloor_MA_mu2base);

Slim_PA_second = ...
    10*log10(Pfloor_PA_mu2red/Pfloor_PA_mu2base);

fprintf('\n===== HIGH-SOURCE-POWER SENSITIVITY LIMITS =====\n');
fprintf('MA first hop  -> %.6f dB\n',Slim_MA_first);
fprintf('MA second hop -> %.6f dB\n',Slim_MA_second);
fprintf('PA first hop  -> %.6f dB\n',Slim_PA_first);
fprintf('PA second hop -> %.6f dB\n',Slim_PA_second);

fprintf('\n===== MAX |ANALYSIS - SIMULATION| =====\n');
fprintf('MA first  = %.4e dB\n', ...
    max(abs(Sana_MA_first-Ssim_MA_first)));
fprintf('MA second = %.4e dB\n', ...
    max(abs(Sana_MA_second-Ssim_MA_second)));
fprintf('PA first  = %.4e dB\n', ...
    max(abs(Sana_PA_first-Ssim_PA_first)));
fprintf('PA second = %.4e dB\n', ...
    max(abs(Sana_PA_second-Ssim_PA_second)));

%% ------------------------------------------------------------------------
% 11) REVISED IEEE-STYLE FIGURE 7
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

%% Zero reference
yline(ax,0,':', ...
    'Color',[0.48 0.48 0.48], ...
    'LineWidth',0.75, ...
    'HandleVisibility','off');

%% Analytical curves
hMA1 = plot(ax, ...
    snrDb,Sana_MA_first,'-', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','MA, first hop');

hMA2 = plot(ax, ...
    snrDb,Sana_MA_second,'--', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','MA, second hop');

hPA1 = plot(ax, ...
    snrDb,Sana_PA_first,':', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','PA, first hop');

hPA2 = plot(ax, ...
    snrDb,Sana_PA_second,'-.', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','PA, second hop');

%% Monte-Carlo markers
% Same marker density and open-circle style as revised Fig. 2.
mk = 1:2:nSNR;                     % every 4 dB

plot(ax, ...
    snrDb(mk),Ssim_MA_first(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(mk),Ssim_MA_second(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(mk),Ssim_PA_first(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(mk),Ssim_PA_second(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

%% Dummy simulation legend handle
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
    'Hop-specific outage change (dB)', ...
    'Interpreter','latex', ...
        'FontSize',10);

xlim(ax,[0 40]);
xticks(ax,0:5:40);

% Automatic y-range, but retain enough room around zero.
allY = [ ...
    Sana_MA_first,Sana_MA_second, ...
    Sana_PA_first,Sana_PA_second, ...
    Ssim_MA_first,Ssim_MA_second, ...
    Ssim_PA_first,Ssim_PA_second];

yMin = min(allY);
yMax = max(allY);
span = max(yMax-yMin,1);

ylim(ax,[min(-0.25,yMin-0.05*span), ...
         max( 0.50,yMax+0.08*span)]);

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
    [hMA1 hMA2 hPA1 hPA2 hSim], ...
    'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',7, ...
    'Box','on');

lgd.LineWidth = 0.55;
lgd.ItemTokenSize = [18 8];

%% Vector export
set(fig,'Renderer','painters');

exportgraphics(fig, ...
    'Fig7_HopSpecific_Fading_Sensitivity_Revised.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'Fig7_HopSpecific_Fading_Sensitivity_Revised.png', ...
    'Resolution',600);

%% ------------------------------------------------------------------------
% 12) APPEND FIGURE-7 RESULTS TO MASTER FILE
% -------------------------------------------------------------------------
save(masterFile, ...
    'M1_base','M2_base','M1_mu1','M2_mu1', ...
    'Pana_MA_first','Pana_MA_second', ...
    'Pana_PA_first','Pana_PA_second', ...
    'Psim_MA_first','Psim_MA_second', ...
    'Psim_PA_first','Psim_PA_second', ...
    'rho_MA_first','rho_MA_second', ...
    'rho_PA_first','rho_PA_second', ...
    'Sana_MA_first','Sana_MA_second', ...
    'Sana_PA_first','Sana_PA_second', ...
    'Ssim_MA_first','Ssim_MA_second', ...
    'Ssim_PA_first','Ssim_PA_second', ...
    'Slim_MA_first','Slim_MA_second', ...
    'Slim_PA_first','Slim_PA_second', ...
    '-append');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function M = normalizationM(alpha,beta,mu,s)
% Unit-mean power normalization:
%
%   M = E[r(U)^2],   U ~ Gamma(mu,1)
%
% where r(.) is the inverse Log-alpha-mu envelope mapping.

r = @(u) ...
    ((1+beta*u/mu).^(1/beta)-s).^(1/alpha) ...
    -(1-s)^(1/alpha);

gpdf = @(u) ...
    u.^(mu-1).*exp(-u)./gamma(mu);

M = integral( ...
    @(u) r(u).^2.*gpdf(u), ...
    0,Inf, ...
    'RelTol',1e-11, ...
    'AbsTol',1e-13);

end

%% ------------------------------------------------------------------------

function [Pbest,rhoBest,Nused] = optimizeRhoAdaptive( ...
    snrDb,Omega1,Omega2,rhoGrid,p,orders,relTol,absTol)

prevBest = NaN;
Pbest = NaN;
rhoBest = NaN;
Nused = orders(end);

for ii = 1:numel(orders)

    N = orders(ii);

    [xGL,wGL] = gaussLegendre(N);

    Pvec = outageGLgrid( ...
        snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL);

    [currBest,idx] = min(Pvec);
    currRho = rhoGrid(idx);

    if ii > 1

        err = abs(currBest-prevBest);
        tol = absTol + relTol*abs(currBest);

        if err <= tol
            Pbest = currBest;
            rhoBest = currRho;
            Nused = N;
            return;
        end
    end

    prevBest = currBest;
    Pbest = currBest;
    rhoBest = currRho;
end

warning(['Gauss-Legendre convergence target was not reached ', ...
         'by N=%d at SNR %.1f dB. Highest-order result is used.'], ...
         orders(end),snrDb);

end

%% ------------------------------------------------------------------------

function Pout = outageGLgrid( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL)
% Finite-order evaluation of Eq. (30) for all rho values.

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

% Transform GL nodes from [-1,1] to [pth,1].
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

Pout = min(max(Pout,realmin),1);

end

%% ------------------------------------------------------------------------

function Pout = directMonteCarloOutage( ...
    Ps,rho,Omega1,Omega2,Z1,Z2,p)
% Direct Eqs. (7)-(10).

X = Omega1*Z1;
Y = Omega2*Z2;

Pin = rho*Ps.*X;
PH = EHmap(Pin,p);

gamma1 = ...
    (1-rho)*Ps.*X/p.sigmaR2;

gamma2 = ...
    PH.*Y/p.sigmaD2;

gammaE2E = ...
    (gamma1.*gamma2)./(gamma1+gamma2+1);

Pout = mean(gammaE2E < p.gammaTh);

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
