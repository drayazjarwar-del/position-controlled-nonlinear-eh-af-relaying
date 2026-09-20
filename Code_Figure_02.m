%% ARCHIVAL FILE: Code_Figure_02.m
% This archival copy preserves the numerical model, parameters, and algorithms
% of the supplied final script. Only dependency/output filenames were standardized
% for the reproducibility package.

%% FIGURE 2: Outage probability vs normalized source transmit power
% Analytical outage from Eq. (30) + direct Monte-Carlo validation
% Nonlinear EH: Eq. (7)
% End-to-end AF SNR: Eq. (10)
% Position optimization: Eqs. (46)-(47)
% Power-splitting optimization: Eqs. (49)-(50)
% High-source-power outage floor: Eq. (34)
% Matched proportional-EH baseline for M6:
%   P_H^lin = eta0*rho*P_S*X, eta0 = Phi'_EH(0) = Pmax*aEH*xi0
%   (no saturation; rho is reoptimized independently at every SNR)
%
% FINAL FIGURE STYLE:
% - Analytical lines pass visibly through open simulation circles.
% - Main analytical lines: 0.90 pt.
% - Asymptotic floors are lighter (0.85 pt) and use a separate color-coded legend.
% - Axis labels use both the physical name and the mathematical symbol.
% - "Normalized source transmit power" is retained consistently as P_S/sigma_R^2.
%
% IMPORTANT:
% The MA coefficients below define ONE explicit reproducible spatial
% realization. If this realization is adopted for the paper, keep it fixed
% for ALL subsequent figures.

clc; close all;

% Resolve this script's own folder so shared data are found reliably.
thisScript = mfilename('fullpath');
if isempty(thisScript)
    thisDir = pwd;
else
    thisDir = fileparts(thisScript);
end
masterFile = fullfile(thisDir,'Master_Numerical_Data.mat');


%% ------------------------------------------------------------------------
% 1) TABLE-II PARAMETERS
% -------------------------------------------------------------------------
c0 = 299792458;
fc = 28e9;
lambda_c = c0/fc;

R0 = 1;                              % bit/s/Hz
gammaTh = 2^(2*R0)-1;                % half-duplex threshold

sigmaR2 = 1;
sigmaD2 = 1;

Pmax = sigmaD2*10^(13/10);           % Pmax/sigma_D^2 = 13 dB
bEH  = sigmaR2*10^(10/10);           % bEH/sigma_R^2 = 10 dB
aEH  = 0.5/sigmaR2;                  % aEH*sigma_R^2 = 0.5
xi0  = 1/(1+exp(aEH*bEH));
eta0 = Pmax*aEH*xi0;                  % Phi'_EH(0), matched linear-EH slope

rhoMin = 0.01;
Nrho = 1001;
rhoGrid = linspace(rhoMin,1-rhoMin,Nrho);

Nmc = 1e6;
snrDb = 0:2:40;

% Default Log-alpha-mu tuple on both hops
alpha1 = 2; beta1 = 1; mu1 = 2; s1 = 0.5;
alpha2 = 2; beta2 = 1; mu2 = 2; s2 = 0.5;

% Numerical convergence settings
relTol = 1e-9;
absTol = 1e-12;
GLorders = [32 64 128 256 512];

%% ------------------------------------------------------------------------
% 2) UNIT-MEAN LOG-alpha-mu NORMALIZATION
% -------------------------------------------------------------------------
r1 = @(u) ((1+beta1*u/mu1).^(1/beta1)-s1).^(1/alpha1) ...
          -(1-s1)^(1/alpha1);

r2 = @(u) ((1+beta2*u/mu2).^(1/beta2)-s2).^(1/alpha2) ...
          -(1-s2)^(1/alpha2);

gpdf1 = @(u) u.^(mu1-1).*exp(-u)./gamma(mu1);
gpdf2 = @(u) u.^(mu2-1).*exp(-u)./gamma(mu2);

M1 = integral(@(u) r1(u).^2.*gpdf1(u),0,Inf, ...
    'RelTol',1e-11,'AbsTol',1e-13);

M2 = integral(@(u) r2(u).^2.*gpdf2(u),0,Inf, ...
    'RelTol',1e-11,'AbsTol',1e-13);

p.alpha1=alpha1; p.beta1=beta1; p.mu1=mu1; p.s1=s1; p.M1=M1;
p.alpha2=alpha2; p.beta2=beta2; p.mu2=mu2; p.s2=s2; p.M2=M2;

p.sigmaR2=sigmaR2;
p.sigmaD2=sigmaD2;

p.Pmax=Pmax;
p.bEH=bEH;
p.aEH=aEH;
p.xi0=xi0;
p.eta0=eta0;
p.gammaTh=gammaTh;

%% ------------------------------------------------------------------------
% 3) EXPLICIT REPRODUCIBLE MA REALIZATION, Eq. (2)
% -------------------------------------------------------------------------
qMA = (0:0.1:4)*lambda_c;
qMA0 = 0;

% Four fixed path weights and direction projections for hop 1
a1 = [-0.75570180-0.42840421j, ...
      -0.19801902-0.04260239j, ...
       0.00920907-0.08588655j, ...
       0.43777523-0.07245455j];

u1 = [ 0.81357359, ...
       0.84949100, ...
      -0.93598546, ...
       0.27135518];

% Four fixed path weights and direction projections for hop 2
a2 = [ 0.03434977+0.33149494j, ...
      -0.04367351-0.56456942j, ...
      -0.55747359+0.40741479j, ...
       0.00173578-0.30252036j];

u2 = [ 0.15892667, ...
      -0.30918135, ...
       0.80069532, ...
      -0.41802376];

% Normalize so that sum_l |a_{i,l}|^2 = 1
a1 = a1/sqrt(sum(abs(a1).^2));
a2 = a2/sqrt(sum(abs(a2).^2));

phase1 = exp(-1j*(2*pi/lambda_c)*(u1(:)*(qMA-qMA0)));
phase2 = exp(-1j*(2*pi/lambda_c)*(u2(:)*(qMA-qMA0)));

psiMA1 = a1(:).' * phase1;
psiMA2 = a2(:).' * phase2;

OmegaMA1 = abs(psiMA1).^2;
OmegaMA2 = abs(psiMA2).^2;

% Normalize fixed reference position q=0 to Omega_i=1
OmegaMA1 = OmegaMA1/OmegaMA1(1);
OmegaMA2 = OmegaMA2/OmegaMA2(1);

[Omega1_MA,idxMA1] = max(OmegaMA1);
[Omega2_MA,idxMA2] = max(OmegaMA2);

qSstar_MA = qMA(idxMA1);
qDstar_MA = qMA(idxMA2);

%% ------------------------------------------------------------------------
% 4) PA REALIZATION, Eqs. (3)-(4)
% -------------------------------------------------------------------------
qPA = 0:0.1:4;                       % m

ellWG = 0.08;                         % dB/m
dPerp1 = 3;
dPerp2 = 3;

qbar1 = 3;
qbar2 = 1;

nPL1 = 2;
nPL2 = 2;

Awg1 = 10.^(-ellWG*qPA/20);
Awg2 = 10.^(-ellWG*qPA/20);

d1 = sqrt(dPerp1^2+(qPA-qbar1).^2);
d2 = sqrt(dPerp2^2+(qPA-qbar2).^2);

% |psi_rad(q)| = 1
OmegaPA1 = Awg1.^2 .* d1.^(-nPL1);
OmegaPA2 = Awg2.^2 .* d2.^(-nPL2);

% Normalize fixed reference position q=0 to Omega_i=1
OmegaPA1 = OmegaPA1/OmegaPA1(1);
OmegaPA2 = OmegaPA2/OmegaPA2(1);

[Omega1_PA,idxPA1] = max(OmegaPA1);
[Omega2_PA,idxPA2] = max(OmegaPA2);

qSstar_PA = qPA(idxPA1);
qDstar_PA = qPA(idxPA2);

%% Fixed-position benchmark
Omega1_FIX = 1;
Omega2_FIX = 1;

fprintf('\n===== OPTIMIZED POSITION GAINS =====\n');

fprintf('MA Omega1*=%.6f (%.3f dB), qS*=%.2f lambda\n', ...
    Omega1_MA,10*log10(Omega1_MA),qSstar_MA/lambda_c);

fprintf('MA Omega2*=%.6f (%.3f dB), qD*=%.2f lambda\n', ...
    Omega2_MA,10*log10(Omega2_MA),qDstar_MA/lambda_c);

fprintf('PA Omega1*=%.6f (%.3f dB), qS*=%.2f m\n', ...
    Omega1_PA,10*log10(Omega1_PA),qSstar_PA);

fprintf('PA Omega2*=%.6f (%.3f dB), qD*=%.2f m\n', ...
    Omega2_PA,10*log10(Omega2_PA),qDstar_PA);

%% ------------------------------------------------------------------------
% 5) ANALYTICAL OUTAGE — EQ. (30)
% -------------------------------------------------------------------------
nSNR = numel(snrDb);

Pana_MA  = zeros(1,nSNR);
Pana_PA  = zeros(1,nSNR);
Pana_FIX = zeros(1,nSNR);

% Matched proportional-EH baselines (optimized MA/PA only)
Pana_LIN_MA = zeros(1,nSNR);
Pana_LIN_PA = zeros(1,nSNR);

rho_MA  = zeros(1,nSNR);
rho_PA  = zeros(1,nSNR);
rho_FIX = zeros(1,nSNR);
rho_LIN_MA = zeros(1,nSNR);
rho_LIN_PA = zeros(1,nSNR);

Nused_MA  = zeros(1,nSNR);
Nused_PA  = zeros(1,nSNR);
Nused_FIX = zeros(1,nSNR);
Nused_LIN_MA = zeros(1,nSNR);
Nused_LIN_PA = zeros(1,nSNR);

fprintf('\n===== ANALYTICAL OUTAGE =====\n');

for k = 1:nSNR

    thisSNR = snrDb(k);

    [Pana_MA(k),rho_MA(k),Nused_MA(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_MA,Omega2_MA, ...
        rhoGrid,p,GLorders,relTol,absTol);

    [Pana_PA(k),rho_PA(k),Nused_PA(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_PA,Omega2_PA, ...
        rhoGrid,p,GLorders,relTol,absTol);

    [Pana_FIX(k),rho_FIX(k),Nused_FIX(k)] = ...
        optimizeRhoAdaptive( ...
        thisSNR,Omega1_FIX,Omega2_FIX, ...
        rhoGrid,p,GLorders,relTol,absTol);

    % Matched proportional-EH baselines. The same optimized positions are
    % used, but rho is reoptimized independently for the linear EH model.
    [Pana_LIN_MA(k),rho_LIN_MA(k),Nused_LIN_MA(k)] = ...
        optimizeRhoAdaptiveLinear( ...
        thisSNR,Omega1_MA,Omega2_MA, ...
        rhoGrid,p,GLorders,relTol,absTol);

    [Pana_LIN_PA(k),rho_LIN_PA(k),Nused_LIN_PA(k)] = ...
        optimizeRhoAdaptiveLinear( ...
        thisSNR,Omega1_PA,Omega2_PA, ...
        rhoGrid,p,GLorders,relTol,absTol);

    fprintf(['SNR=%2d dB | NL-MA %.4e rho=%.3f | ', ...
             'NL-PA %.4e rho=%.3f | Fixed %.4e rho=%.3f | ', ...
             'LIN-MA %.4e rho=%.3f | LIN-PA %.4e rho=%.3f\n'], ...
             thisSNR, ...
             Pana_MA(k),rho_MA(k), ...
             Pana_PA(k),rho_PA(k), ...
             Pana_FIX(k),rho_FIX(k), ...
             Pana_LIN_MA(k),rho_LIN_MA(k), ...
             Pana_LIN_PA(k),rho_LIN_PA(k));
end

%% ------------------------------------------------------------------------
% 6) MONTE-CARLO VALIDATION — DIRECT EQS. (7) AND (10)
% -------------------------------------------------------------------------
rng(20260910,'twister');

U1mc = randg(mu1,Nmc,1);
U2mc = randg(mu2,Nmc,1);

Z1mc = r1(U1mc).^2/M1;
Z2mc = r2(U2mc).^2/M2;

Psim_MA  = zeros(1,nSNR);
Psim_PA  = zeros(1,nSNR);
Psim_FIX = zeros(1,nSNR);

for k = 1:nSNR

    Ps = sigmaR2*10^(snrDb(k)/10);

    Psim_MA(k) = directMonteCarloOutage( ...
        Ps,rho_MA(k),Omega1_MA,Omega2_MA,Z1mc,Z2mc,p);

    Psim_PA(k) = directMonteCarloOutage( ...
        Ps,rho_PA(k),Omega1_PA,Omega2_PA,Z1mc,Z2mc,p);

    Psim_FIX(k) = directMonteCarloOutage( ...
        Ps,rho_FIX(k),Omega1_FIX,Omega2_FIX,Z1mc,Z2mc,p);
end

%% ------------------------------------------------------------------------
% 7) ASYMPTOTIC OUTAGE FLOORS — EQ. (34)
% -------------------------------------------------------------------------
Pfloor_MA = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_MA), ...
    alpha2,beta2,mu2,s2,M2);

Pfloor_PA = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_PA), ...
    alpha2,beta2,mu2,s2,M2);

Pfloor_FIX = FZ( ...
    gammaTh*sigmaD2/(Pmax*Omega2_FIX), ...
    alpha2,beta2,mu2,s2,M2);

fprintf('\n===== ASYMPTOTIC FLOORS =====\n');
fprintf('MA    = %.6e\n',Pfloor_MA);
fprintf('PA    = %.6e\n',Pfloor_PA);
fprintf('Fixed = %.6e\n',Pfloor_FIX);

fprintf('\n===== MAX ANALYSIS-SIMULATION DIFFERENCE =====\n');
fprintf('MA    = %.4e\n',max(abs(Pana_MA-Psim_MA)));
fprintf('PA    = %.4e\n',max(abs(Pana_PA-Psim_PA)));
fprintf('Fixed = %.4e\n',max(abs(Pana_FIX-Psim_FIX)));

%% ------------------------------------------------------------------------
% 8) REVISED IEEE-STYLE FIGURE 2
% -------------------------------------------------------------------------
% Restrained colors + line-style differences for grayscale robustness
cMA  = [1.00 0.00 0.00];
cPA  = [0.00 0.00 1.00];
cFIX = [0.00 0.00 0.00];
cSIM = [0.00 0.00 0.00];

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 3.50 2.55]);

ax = axes(fig);
hold(ax,'on');
box(ax,'on');

%% Analytical curves: reduced line width
hMA = semilogy(ax, ...
    snrDb,Pana_MA,'-', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','Optimized MA');

hPA = semilogy(ax, ...
    snrDb,Pana_PA,'--', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','Optimized PA');

hFIX = semilogy(ax, ...
    snrDb,Pana_FIX,':', ...
    'Color',cFIX, ...
    'LineWidth',0.90, ...
    'DisplayName','Fixed positions');

%% Matched proportional-EH baselines (analytical; no MC validation markers)
% Keep architecture colors, but use styles distinct from both the nonlinear
% curves and the dash-dot asymptotic floors. Sparse +/x markers identify
% these analytical baselines and remain distinct from open-circle MC markers.
linMk = 1:3:nSNR;                     % every 6 dB

hLinMA = semilogy(ax, ...
    snrDb,Pana_LIN_MA,'--', ...
    'Color',cMA, ...
    'LineWidth',0.85, ...
    'Marker','+', ...
    'MarkerIndices',linMk, ...
    'MarkerSize',3.2, ...
    'DisplayName','Linear EH: MA');

hLinPA = semilogy(ax, ...
    snrDb,Pana_LIN_PA,':', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'Marker','x', ...
    'MarkerIndices',linMk, ...
    'MarkerSize',3.2, ...
    'DisplayName','Linear EH: PA');

%% Asymptotic floors: lighter than main curves
semilogy(ax, ...
    snrDb,Pfloor_MA*ones(size(snrDb)),'-.', ...
    'Color',cMA, ...
    'LineWidth',0.85, ...
    'HandleVisibility','off');

semilogy(ax, ...
    snrDb,Pfloor_PA*ones(size(snrDb)),'-.', ...
    'Color',cPA, ...
    'LineWidth',0.85, ...
    'HandleVisibility','off');

semilogy(ax, ...
    snrDb,Pfloor_FIX*ones(size(snrDb)),'-.', ...
    'Color',cFIX, ...
    'LineWidth',0.85, ...
    'HandleVisibility','off');

%% Simulation markers
% Open black circles are reserved exclusively for Monte-Carlo simulation.
% The analytical lines remain visible through the unfilled markers.
mk = 1:2:nSNR;                        % every 4 dB

semilogy(ax, ...
    snrDb(mk),Psim_MA(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

semilogy(ax, ...
    snrDb(mk),Psim_PA(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

semilogy(ax, ...
    snrDb(mk),Psim_FIX(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

%% Dummy simulation handle for the main legend
hSim = semilogy(ax, ...
    NaN,NaN,'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'DisplayName','Simulation');

%% Axis labels: physical name + mathematical symbol
xlabel(ax, ...
    'Normalized source transmit power (dB)', ...
    'Interpreter','latex', ...
    'FontSize',10);

ylabel(ax, ...
    'Outage probability', ...
    'Interpreter','latex', ...
    'FontSize',10);

xlim(ax,[0 40]);
ylim(ax,[1e-3 1]);

xticks(ax,0:5:40);
yticks(ax,[1e-3 1e-2 1e-1 1]);

%% Clean IEEE-style axes/grid
set(ax, ...
    'YScale','log', ...
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

%% Main legend: nonlinear curves + proportional-EH baselines + MC
% Use two columns to keep the legend compact and place it in the upper-right
% region, away from the high-power linear-EH separation that motivates M6.
lgdMain = legend(ax, ...
    [hMA hPA hFIX hLinMA hLinPA hSim], ...
    {'Optimized MA','Optimized PA','Fixed positions', ...
     'Linear EH: MA','Linear EH: PA','Simulation'}, ...
    'Location','northeast', ...
    'NumColumns',2, ...
    'Interpreter','latex', ...
    'FontSize',6.2, ...
    'Box','on');

lgdMain.LineWidth = 0.55;
lgdMain.ItemTokenSize = [15 8];
lgdMain.AutoUpdate = 'off';

%% Second legend: architecture-specific asymptotic floors
% Place a separate vertical legend-style box on the LEFT side of the axes,
% in the open region between the upper PA/fixed floors and the lower MA
% floor.  An annotation-based box is used because MATLAB normally allows
% only one ordinary legend per axes and this method is reliable in vector
% export across MATLAB releases.
drawnow;
ax.Units = 'normalized';
pAx = ax.Position;

% Position relative to the plotting axes so that it scales with the figure.
% [left bottom width height] in normalized figure coordinates.
floorW = 0.33*pAx(3);
floorH = 0.23*pAx(4);
floorX = pAx(1) + 0.045*pAx(3);
floorY = pAx(2) + 0.210*pAx(4);
pFloor = [floorX floorY floorW floorH];

% White legend background and 0.55-pt border.
annotation(fig,'rectangle',pFloor, ...
    'Color',[0 0 0], ...
    'FaceColor',[1 1 1], ...
    'LineWidth',0.55);

% Three vertically stacked floor entries.
padX   = 0.055*floorW;
tokenW = 0.27*floorW;
gapTX  = 0.040*floorW;
rowH   = floorH/3;

floorLabels = {'MA floor','PA floor','Fixed floor'};
floorColors = [cMA; cPA; cFIX];

for ii = 1:3
    % Top-to-bottom row order: MA, PA, Fixed.
    yMid = floorY + floorH - (ii-0.5)*rowH;
    x0   = floorX + padX;

    annotation(fig,'line', ...
        [x0 x0+tokenW], ...
        [yMid yMid], ...
        'Color',floorColors(ii,:), ...
        'LineStyle','-.', ...
        'LineWidth',0.85);

    annotation(fig,'textbox', ...
        [x0+tokenW+gapTX, yMid-0.42*rowH, ...
         floorW-padX-tokenW-gapTX-0.02*floorW, 0.84*rowH], ...
        'String',floorLabels{ii}, ...
        'Interpreter','latex', ...
        'FontName','Times New Roman', ...
        'FontSize',7, ...
        'EdgeColor','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'Margin',0);
end

% Restore the plotting axes as current axes.
axes(ax);

%% Vector renderer and export
set(fig,'Renderer','painters');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_02.pdf'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(thisDir,'Figure_02.png'), ...
    'Resolution',600);

%% ------------------------------------------------------------------------
% 9) SAVE MASTER NUMERICAL DATA FOR LATER FIGURES
% -------------------------------------------------------------------------
save(masterFile, ...
    'snrDb', ...
    'Pana_MA','Pana_PA','Pana_FIX', ...
    'Pana_LIN_MA','Pana_LIN_PA', ...
    'Psim_MA','Psim_PA','Psim_FIX', ...
    'rho_MA','rho_PA','rho_FIX', ...
    'rho_LIN_MA','rho_LIN_PA', ...
    'Pfloor_MA','Pfloor_PA','Pfloor_FIX', ...
    'OmegaMA1','OmegaMA2', ...
    'OmegaPA1','OmegaPA2', ...
    'Omega1_MA','Omega2_MA', ...
    'Omega1_PA','Omega2_PA', ...
    'qMA','qPA', ...
    'qSstar_MA','qDstar_MA', ...
    'qSstar_PA','qDstar_PA', ...
    'a1','a2','u1','u2', ...
    'M1','M2','eta0');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

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
        tol = absTol+relTol*abs(currBest);

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

warning( ...
    ['Gauss-Legendre convergence target was not reached ', ...
     'by N=%d at SNR %.1f dB. Highest-order result is used.'], ...
     orders(end),snrDb);

end

%% ------------------------------------------------------------------------

function Pout = outageGLgrid( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL)
% Finite-order evaluation of Eq. (30) for all rho-grid values.

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

% Map Gauss-Legendre nodes from [-1,1] to [pth,1].
v = (1+pa)/2 + ((1-pa)/2).*xGL(:).';

% Avoid numerical inverse-gamma evaluation at exactly 0 or 1.
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

Pout = min(max(Pout,0),1);

end

%% ------------------------------------------------------------------------

function [Pbest,rhoBest,Nused] = optimizeRhoAdaptiveLinear( ...
    snrDb,Omega1,Omega2,rhoGrid,p,orders,relTol,absTol)
% Adaptive Gauss-Legendre evaluation for the matched proportional-EH
% baseline P_H^lin = eta0*rho*P_S*X. rho is optimized independently
% from the nonlinear-EH design.

prevBest = NaN;
Pbest = NaN;
rhoBest = NaN;
Nused = orders(end);

for ii = 1:numel(orders)

    N = orders(ii);
    [xGL,wGL] = gaussLegendre(N);

    Pvec = outageGLgridLinear( ...
        snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL);

    [currBest,idx] = min(Pvec);
    currRho = rhoGrid(idx);

    if ii > 1

        err = abs(currBest-prevBest);
        tol = absTol+relTol*abs(currBest);

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

warning( ...
    ['Linear-EH Gauss-Legendre convergence target was not reached ', ...
     'by N=%d at SNR %.1f dB. Highest-order result is used.'], ...
     orders(end),snrDb);

end

%% ------------------------------------------------------------------------

function Pout = outageGLgridLinear( ...
    snrDb,rhoGrid,Omega1,Omega2,p,xGL,wGL)
% Matched proportional-EH counterpart of Eq. (30):
%   P_H^lin = eta0*rho*P_S*X,
% where eta0 = Phi'_EH(0) = Pmax*aEH*xi0. No saturation clipping.

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

% Map Gauss-Legendre nodes from [-1,1] to [pth,1].
v = (1+pa)/2 + ((1-pa)/2).*xGL(:).';
v = min(v,1-eps);
v = max(v,realmin);

u = gammaincinv(v,p.mu1,'lower');

r = ((1+p.beta1*u/p.mu1).^(1/p.beta1)-p.s1) ...
    .^(1/p.alpha1) ...
    -(1-p.s1)^(1/p.alpha1);

x = Omega1*r.^2/p.M1;

% Matched proportional EH model; deliberately unbounded.
PH = p.eta0*rho(active).*Ps.*x;
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

function Pout = directMonteCarloOutage( ...
    Ps,rho,Omega1,Omega2,Z1,Z2,p)
% Direct system simulation from Eqs. (7)-(10).

X = Omega1*Z1;
Y = Omega2*Z2;

Pin = rho*Ps.*X;
PH = EHmap(Pin,p);

gamma1 = (1-rho)*Ps.*X/p.sigmaR2;
gamma2 = PH.*Y/p.sigmaD2;

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
% Eq. (13), using the regularized lower incomplete gamma.

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
