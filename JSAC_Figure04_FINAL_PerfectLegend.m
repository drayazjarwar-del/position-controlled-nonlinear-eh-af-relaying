%% FIGURE 4: EH-induced source-to-destination positioning transition
% Analytical result + direct Monte-Carlo validation
%
% The plotted metric is
%
%   D^(chi) = 10 log10(Pout,dest-only / Pout,source-only)
%
% which is equivalently the difference (in dB) between the outage reduction
% produced by source-only positioning and that produced by destination-only
% positioning, both relative to the same fixed-position benchmark.
%
% Therefore:
%   D^(chi) > 0  -> source-side positioning gives the larger reliability gain
%   D^(chi) < 0  -> destination-side positioning gives the larger gain
%
% IMPORTANT:
% 1) Run JSAC_Figure02_FINAL.m first.
% 2) This script loads JSAC_Master_Numerical_Data.mat, so Fig. 4 uses the
%    SAME MA/PA spatial realization as Figs. 2 and 3.
% 3) The power-splitting ratio is re-optimized separately for source-only
%    and destination-only positioning, as stated in the manuscript.
% 4) The noncircular markers identify the normalized-source-transmit-power points at which
%    Pr{P_H >= 0.9 Pmax} reaches 0.1, 0.5, and 0.9 for the FULLY optimized
%    MA/PA design of Fig. 2. This makes the saturation markers unambiguous.
%
% IEEE-style formatting:
% - analytical lines: 0.90 pt
% - simulation circles: open/transparent
% - zero-reference line: thin gray dotted line
% - saturation markers: square / diamond / triangle
% - one-column width: 3.50 in

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
    'rho_MA','rho_PA', ...
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

Nmc = 1e6;

% Default Log-alpha-mu tuple
alpha1 = 2; beta1 = 1; mu1 = 2; s1 = 0.5;
alpha2 = 2; beta2 = 1; mu2 = 2; s2 = 0.5;

% Adaptive Gauss-Legendre settings for Eq. (30)
relTol = 1e-9;
absTol = 1e-12;
GLorders = [32 64 128 256 512];

p.alpha1=alpha1; p.beta1=beta1; p.mu1=mu1; p.s1=s1; p.M1=M1;
p.alpha2=alpha2; p.beta2=beta2; p.mu2=mu2; p.s2=s2; p.M2=M2;
p.sigmaR2=sigmaR2; p.sigmaD2=sigmaD2;
p.Pmax=Pmax; p.bEH=bEH; p.aEH=aEH; p.xi0=xi0;
p.gammaTh=gammaTh;

Omega1_FIX = 1;
Omega2_FIX = 1;

%% ------------------------------------------------------------------------
% 3) ANALYTICAL SOURCE-ONLY / DESTINATION-ONLY OUTAGE
% -------------------------------------------------------------------------
nSNR = numel(snrDb);

% Source-only positioning:
%   source uses Omega1*, destination remains fixed at Omega2 = 1.
Psrc_MA  = zeros(1,nSNR);
Psrc_PA  = zeros(1,nSNR);
rhoSrc_MA = zeros(1,nSNR);
rhoSrc_PA = zeros(1,nSNR);

% Destination-only positioning:
%   source remains fixed at Omega1 = 1, destination uses Omega2*.
Pdst_MA  = zeros(1,nSNR);
Pdst_PA  = zeros(1,nSNR);
rhoDst_MA = zeros(1,nSNR);
rhoDst_PA = zeros(1,nSNR);

fprintf('\n===== FIG. 4: ANALYTICAL SOURCE-/DESTINATION-ONLY OUTAGE =====\n');

for k = 1:nSNR

    thisSNR = snrDb(k);

    [Psrc_MA(k),rhoSrc_MA(k)] = ...
        optimizeRhoAdaptive(thisSNR,Omega1_MA,Omega2_FIX, ...
                            rhoGrid,p,GLorders,relTol,absTol);

    [Pdst_MA(k),rhoDst_MA(k)] = ...
        optimizeRhoAdaptive(thisSNR,Omega1_FIX,Omega2_MA, ...
                            rhoGrid,p,GLorders,relTol,absTol);

    [Psrc_PA(k),rhoSrc_PA(k)] = ...
        optimizeRhoAdaptive(thisSNR,Omega1_PA,Omega2_FIX, ...
                            rhoGrid,p,GLorders,relTol,absTol);

    [Pdst_PA(k),rhoDst_PA(k)] = ...
        optimizeRhoAdaptive(thisSNR,Omega1_FIX,Omega2_PA, ...
                            rhoGrid,p,GLorders,relTol,absTol);

    fprintf(['SNR=%2d dB | MA: Psrc %.4e, Pdst %.4e | ', ...
             'PA: Psrc %.4e, Pdst %.4e\n'], ...
             thisSNR,Psrc_MA(k),Pdst_MA(k), ...
             Psrc_PA(k),Pdst_PA(k));
end

%% ------------------------------------------------------------------------
% 4) RELATIVE POSITIONING-GAIN METRIC
% -------------------------------------------------------------------------
% Difference between the source-only and destination-only outage reductions:
%
% [10log10(Pfix/Psrc)] - [10log10(Pfix/Pdst)]
%       = 10log10(Pdst/Psrc)
%
Dana_MA = 10*log10(Pdst_MA./Psrc_MA);
Dana_PA = 10*log10(Pdst_PA./Psrc_PA);

%% ------------------------------------------------------------------------
% 5) MONTE-CARLO VALIDATION
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

PsrcSim_MA = zeros(1,nSNR);
PdstSim_MA = zeros(1,nSNR);
PsrcSim_PA = zeros(1,nSNR);
PdstSim_PA = zeros(1,nSNR);

for k = 1:nSNR

    Ps = sigmaR2*10^(snrDb(k)/10);

    PsrcSim_MA(k) = directMonteCarloOutage( ...
        Ps,rhoSrc_MA(k),Omega1_MA,Omega2_FIX,Z1mc,Z2mc,p);

    PdstSim_MA(k) = directMonteCarloOutage( ...
        Ps,rhoDst_MA(k),Omega1_FIX,Omega2_MA,Z1mc,Z2mc,p);

    PsrcSim_PA(k) = directMonteCarloOutage( ...
        Ps,rhoSrc_PA(k),Omega1_PA,Omega2_FIX,Z1mc,Z2mc,p);

    PdstSim_PA(k) = directMonteCarloOutage( ...
        Ps,rhoDst_PA(k),Omega1_FIX,Omega2_PA,Z1mc,Z2mc,p);
end

Dsim_MA = 10*log10(PdstSim_MA./PsrcSim_MA);
Dsim_PA = 10*log10(PdstSim_PA./PsrcSim_PA);

%% ------------------------------------------------------------------------
% 6) ZERO-CROSSING / TRANSITION SNR
% -------------------------------------------------------------------------
snrCross_MA = firstZeroCrossing(snrDb,Dana_MA);
snrCross_PA = firstZeroCrossing(snrDb,Dana_PA);

fprintf('\n===== SOURCE-TO-DESTINATION TRANSITION =====\n');

if isnan(snrCross_MA)
    fprintf('MA: no zero crossing in the simulated SNR range.\n');
else
    fprintf('MA zero crossing = %.3f dB\n',snrCross_MA);
end

if isnan(snrCross_PA)
    fprintf('PA: no zero crossing in the simulated SNR range.\n');
else
    fprintf('PA zero crossing = %.3f dB\n',snrCross_PA);
end

%% ------------------------------------------------------------------------
% 7) EH-SATURATION MARKERS
% -------------------------------------------------------------------------
% We use the fully optimized designs from Fig. 2:
%   (Omega1_MA, Omega2_MA, rho_MA)
%   (Omega1_PA, Omega2_PA, rho_PA)
%
% Only Omega1 and rho enter Pr{P_H >= 0.9 Pmax}, but rho is the joint
% outage-optimal rho from Fig. 2.

tauSat = 0.9;
satTargets = [0.1 0.5 0.9];

pSat_MA = zeros(1,nSNR);
pSat_PA = zeros(1,nSNR);

for k = 1:nSNR

    pSat_MA(k) = saturationProbability( ...
        snrDb(k),rho_MA(k),Omega1_MA,tauSat,p);

    pSat_PA(k) = saturationProbability( ...
        snrDb(k),rho_PA(k),Omega1_PA,tauSat,p);
end

satSNR_MA = nan(size(satTargets));
satSNR_PA = nan(size(satTargets));

satD_MA = nan(size(satTargets));
satD_PA = nan(size(satTargets));

for j = 1:numel(satTargets)

    satSNR_MA(j) = thresholdCrossing( ...
        snrDb,pSat_MA,satTargets(j));

    satSNR_PA(j) = thresholdCrossing( ...
        snrDb,pSat_PA,satTargets(j));

    if ~isnan(satSNR_MA(j))
        satD_MA(j) = interp1(snrDb,Dana_MA, ...
            satSNR_MA(j),'pchip');
    end

    if ~isnan(satSNR_PA(j))
        satD_PA(j) = interp1(snrDb,Dana_PA, ...
            satSNR_PA(j),'pchip');
    end
end

fprintf('\n===== Pr{P_H >= 0.9 Pmax} SATURATION MARKERS =====\n');

for j = 1:numel(satTargets)
    fprintf('Target %.1f | MA %.3f dB | PA %.3f dB\n', ...
        satTargets(j),satSNR_MA(j),satSNR_PA(j));
end

fprintf('\n===== MAX |ANALYSIS - SIMULATION| FOR D =====\n');
fprintf('MA = %.4e dB\n',max(abs(Dana_MA-Dsim_MA)));
fprintf('PA = %.4e dB\n',max(abs(Dana_PA-Dsim_PA)));

%% ------------------------------------------------------------------------
% 8) REVISED IEEE-STYLE FIGURE 4
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
hZero = yline(ax,0,':', ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',0.85, ...
    'HandleVisibility','off');

%% Analytical curves
hMA = plot(ax, ...
    snrDb,Dana_MA,'-', ...
    'Color',cMA, ...
    'LineWidth',0.90, ...
    'DisplayName','MA');

hPA = plot(ax, ...
    snrDb,Dana_PA,'--', ...
    'Color',cPA, ...
    'LineWidth',0.90, ...
    'DisplayName','PA');

%% Monte-Carlo circles
mk = 1:2:nSNR;                       % every 4 dB

plot(ax, ...
    snrDb(mk),Dsim_MA(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

plot(ax, ...
    snrDb(mk),Dsim_PA(mk),'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'HandleVisibility','off');

%% Saturation markers
% square:  Pr = 0.1
% diamond: Pr = 0.5
% triangle:Pr = 0.9
satMarkers = {'s','d','^'};

for j = 1:numel(satTargets)

    if ~isnan(satSNR_MA(j))
        plot(ax,satSNR_MA(j),satD_MA(j),satMarkers{j}, ...
            'LineStyle','none', ...
            'MarkerSize',5.0, ...
            'MarkerFaceColor','none', ...
            'MarkerEdgeColor',cSIM, ...
            'LineWidth',1.00, ...
            'HandleVisibility','off');
    end

    if ~isnan(satSNR_PA(j))
        plot(ax,satSNR_PA(j),satD_PA(j),satMarkers{j}, ...
            'LineStyle','none', ...
            'MarkerSize',5.0, ...
            'MarkerFaceColor','none', ...
            'MarkerEdgeColor',cSIM, ...
            'LineWidth',1.00, ...
            'HandleVisibility','off');
    end
end

%% Dummy simulation handle
hSim = plot(ax, ...
    NaN,NaN,'o', ...
    'LineStyle','none', ...
    'MarkerSize',3.8, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'LineWidth',0.70, ...
    'DisplayName','Simulation');

%% Y-axis range
yl = [min([Dana_MA,Dana_PA,Dsim_MA,Dsim_PA])-0.25, ...
      max([Dana_MA,Dana_PA,Dsim_MA,Dsim_PA])+0.25];

% Keep a sensible visual range if data are modest.
yl(1) = min(yl(1),-0.5);
yl(2) = max(yl(2), 0.5);

ylim(ax,yl);

%% Axis labels
xlabel(ax, ...
    'Normalized source transmit power (dB)', ...
    'Interpreter','latex', ...
        'FontSize',10);

ylabel(ax, ...
    'Relative positioning gain (dB)', ...
    'Interpreter','latex', ...
        'FontSize',10);

xlim(ax,[0 40]);
xticks(ax,0:5:40);

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

%% Legends
lgd = legend(ax,[hMA hPA hSim], ...
    'Location','northeast', ...
    'Interpreter','latex', ...
    'FontSize',7, ...
    'Box','on');

lgd.LineWidth = 0.55;
lgd.ItemTokenSize = [20 8];

% Second compact legend for the EH-saturation probability milestones:
% square / diamond / triangle correspond to
% Pr{P_H >= 0.9 P_max} = 0.1 / 0.5 / 0.9, respectively.
axSat = axes(fig, ...
    'Position',ax.Position, ...
    'Color','none', ...
    'XColor','none', ...
    'YColor','none', ...
    'Visible','off', ...
    'HitTest','off');

hold(axSat,'on');

hSat1 = plot(axSat,nan,nan,'s', ...
    'LineStyle','none', ...
    'Color',cSIM, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'MarkerSize',5.0, ...
    'LineWidth',0.90);

hSat2 = plot(axSat,nan,nan,'d', ...
    'LineStyle','none', ...
    'Color',cSIM, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'MarkerSize',5.0, ...
    'LineWidth',0.90);

hSat3 = plot(axSat,nan,nan,'^', ...
    'LineStyle','none', ...
    'Color',cSIM, ...
    'MarkerFaceColor','none', ...
    'MarkerEdgeColor',cSIM, ...
    'MarkerSize',5.0, ...
    'LineWidth',0.90);

lgdSat = legend(axSat,[hSat1 hSat2 hSat3], ...
    {'$\Pr\{P_H \geq 0.9P_{\max}\}=0.1$', ...
     '$\Pr\{P_H \geq 0.9P_{\max}\}=0.5$', ...
     '$\Pr\{P_H \geq 0.9P_{\max}\}=0.9$'}, ...
    'Interpreter','latex', ...
    'FontSize',6.7, ...
    'Box','on');

lgdSat.LineWidth = 0.55;
lgdSat.ItemTokenSize = [9 8];
lgdSat.Units = 'normalized';

% Compact lower-left placement to avoid overlap with the MA curve and
% to remove the excessive blank space inside the second legend.
satLegendX = ax.Position(1) + 0.012*ax.Position(3);
satLegendY = ax.Position(2) + 0.015*ax.Position(4);
satLegendW = 0.43*ax.Position(3);
satLegendH = 0.215*ax.Position(4);
lgdSat.Position = [satLegendX, satLegendY, satLegendW, satLegendH];

%% Export
set(fig,'Renderer','painters');

exportgraphics(fig, ...
    'Fig4_Positioning_Transition_Analysis_Simulation_Revised.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'Fig4_Positioning_Transition_Analysis_Simulation_Revised.png', ...
    'Resolution',600);

%% ------------------------------------------------------------------------
% 9) APPEND FIGURE-4 DATA TO MASTER FILE
% -------------------------------------------------------------------------
save(masterFile, ...
    'Psrc_MA','Pdst_MA','Psrc_PA','Pdst_PA', ...
    'PsrcSim_MA','PdstSim_MA','PsrcSim_PA','PdstSim_PA', ...
    'rhoSrc_MA','rhoDst_MA','rhoSrc_PA','rhoDst_PA', ...
    'Dana_MA','Dana_PA','Dsim_MA','Dsim_PA', ...
    'snrCross_MA','snrCross_PA', ...
    'pSat_MA','pSat_PA','satTargets', ...
    'satSNR_MA','satSNR_PA','satD_MA','satD_PA', ...
    '-append');

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
% Finite-order evaluation of Eq. (30) over the entire rho grid.

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

% Map standard GL nodes from [-1,1] to [pth,1].
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

Pout = min(max(Pout,0),1);

end

%% ------------------------------------------------------------------------

function Pout = directMonteCarloOutage( ...
    Ps,rho,Omega1,Omega2,Z1,Z2,p)
% Direct evaluation of Eqs. (7)-(10).

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

function prob = saturationProbability(snrDb,rho,Omega1,tau,p)
% Exact probability Pr{P_H >= tau Pmax}.
%
% Invert the sigmoidal EH function:
%   PH/Pmax = tau
% to obtain the required RF-input threshold.

Ps = p.sigmaR2*10^(snrDb/10);

L = p.xi0 + (1-p.xi0)*tau;

PinTau = ...
    p.bEH + (1/p.aEH)*log(L/(1-L));

zth = PinTau/(rho*Ps*Omega1);

prob = 1 - FZ( ...
    zth,p.alpha1,p.beta1,p.mu1,p.s1,p.M1);

prob = min(max(prob,0),1);

end

%% ------------------------------------------------------------------------

function xcross = firstZeroCrossing(x,y)
% First zero crossing, linearly interpolated.

xcross = NaN;

for k = 1:numel(x)-1

    if y(k) == 0
        xcross = x(k);
        return;
    end

    if y(k)*y(k+1) < 0

        xcross = ...
            x(k) - y(k)*(x(k+1)-x(k))/(y(k+1)-y(k));

        return;
    end
end

if y(end) == 0
    xcross = x(end);
end

end

%% ------------------------------------------------------------------------

function xcross = thresholdCrossing(x,y,target)
% First x where a monotone/nondecreasing trajectory reaches target.
% Linear interpolation is used between the two bracketing samples.

xcross = NaN;

for k = 1:numel(x)-1

    y1 = y(k)-target;
    y2 = y(k+1)-target;

    if y1 == 0
        xcross = x(k);
        return;
    end

    if y1*y2 < 0 || y2 == 0

        xcross = ...
            x(k) + (target-y(k)) * ...
            (x(k+1)-x(k))/(y(k+1)-y(k));

        return;
    end
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
% Eq. (13) in terms of the regularized lower incomplete gamma.

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
