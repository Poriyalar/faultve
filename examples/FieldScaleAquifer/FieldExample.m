%% Malay Basin Field Scale Simulation Example
% The J area, a conceptual CO2 storage site within the Malay Basin 
% offshore Peninsular Malaysia,is used to illustrate the fault leakage 
% function coupled with VE models for field-scale fault leakage quantification.
% Top of the reservoir lies at depths between 1500 m and 2500 m below 
% a seafloor depthof 70 m. The reservoir pressure is assumed hydrostatic 
% with a temperature profile following a geothermal gradient of 50°C/km
% and a seafloor temperature of 24°C.A 26km-long fault is introduced 
% into this grid and the fault properties are provided. The injection well 
% islocated 12km away from the fault.

% This example is described in more detail in Section 3.4 of the manuscript
% "Rapid Fault Leakage Modeling for CO2 Storage in Saline Aquifers".
% DOI - https://doi.org/10.31223/X5S12N
%% Clear workspace and figures
clear; clc; close all;

%% add modules
mrstModule add faultve
mrstModule add coarsegrid
mrstModule add ad-core
mrstModule add ad-props
gravity on;

%% Load and process data
% the grid (grdecl) file is present in this folder. This file is used to
% generate the 3D grid and VE grid subsequently.
grdecl              = readGRDECL('JENERA.GRDECL');
G                   = processGRDECL(grdecl);
G                   = computeGeometry(G);
[Gt, ~, transMult]  = topSurfaceGrid(G);
% The reservoir rock is created with a 50 [md] permeability and 0.2 [-] porosity
rock                = makeRock(G,50 * milli * darcy, 0.2);
rock2D              = averageRock(rock, Gt);

%% creating fluid and rock fluid parameters (reservoir)
g       = gravity;
rhow    = 1001;                                 % water density [kg m^3]
co2     = CO2props();                           % CO2 property functions
p_ref   = 30 *mega;                             % reference pressure in [Pa]
t_ref   = 94+273.15;                            % reference temperature in [K]
co2_rho = 389.7;                                % CO2 density in [kg m^-3]
co2_c   = 0.0;                                  % CO2 compressibility in [Pa^-1]
wat_c   = 0.0;                                  % water compressibility in [Pa^-1]
c_rock  = 4.35e-5 / barsa;                      % rock compressibility in [Pa^-1]
srw     = 0.3;                                  % residual water [-]
src     = 0.2;                                  % residual CO2 [-]
pe      = 0.0;                                  % capillary entry pressure in [Pa]
muw     = 3.13e-4;                              % brine viscosity in [Pa s]
muco2   = 3.21e-5;                              % CO2 viscosity in [Pa s]
% The PVT behavior of the injected CO2 is assumed to be given by an
% equation state, whereas the brine is incompressible.
% Please consult the documentation for the makeVEFluid routine for a
% description of parameters you can use to set up the various fluid models
% implemented in MRST-co2lab.
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^1; % uses CO2 saturation
fluid   = makeVEFluid(Gt, rock, 'sharp interface'       , ...
               'co2_mu_ref'  , muco2                    , ...
               'wat_mu_ref'  , muw                      , ...
               'co2_rho_ref' , co2_rho                  , ...
               'wat_rho_ref' , rhow                     , ...
               'co2_rho_pvt' , [co2_c, p_ref]           , ...
               'wat_rho_pvt' , [wat_c, p_ref]           , ...
               'residual'    , [srw, src]               , ...
               'pvMult_p_ref', p_ref                    , ...
               'pvMult_fac'  , c_rock                   , ...
               'invPc3D'     , invPc3D                  , ...
               'kr3D'        , kr3D                     , ...
               'transMult'   , transMult);
%% creating the fault block
% The gridblocks  and faces containing the fault block are identified 
% fblock contains the faces that need to be blocked due to fault core
% cblock contains the cells that hold the fault leakage terms
fblock = addSealingFacesX(Gt,'x_range',[425000, 425000],'y_range',[717000, 733000]);
cblock = find(Gt.cells.centroids(1:Gt.cells.num,1) >= 424800 & ...
    Gt.cells.centroids(1:Gt.cells.num,1) <= 425000 & ...
    Gt.cells.centroids(1:Gt.cells.num,2) >= 717000 & ...
    Gt.cells.centroids(1:Gt.cells.num,2) < 733000);

%% creating well, initial state, boundary state and schedule
% The simulation will consist of two periods: during the first 30 years,
% CO2 is injected as constant rate from the single injector. This period is
% simulated with a time step of 1 year. We then simulate a post-injection
% period of 970 years using time steps of 10 years.

% initial conditions
initState.pressure  = rhow * g(3) * Gt.cells.z;
initState.s         = repmat([1, 0], Gt.cells.num, 1);
initState.sGmax     = initState.s(:,2);
inj                 = 1.0;              % injection rate is 1 Mt/year   
inj_yr              = inj*1000000*1000; % convert to kg/yr
inj_yr_s            = inj_yr/year;      % convert to kg/sec
inj_rate            = inj_yr_s/co2_rho; % convert to m3/s

% Well information
wellCell    = findEnclosingCell(Gt, [410100,721700]);
W           = [];
W           = addWellVE(W, Gt, rock2D, wellCell, 'name', 'injector', ...
                'type', 'rate', ...  % inject at constant rate
                'val', inj_rate, ... % volumetric injection rate
                'comp_i', [0 1]);    % inject CO2, not water)

W2D         = W;

% Find bounding indices in 2D grid
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
j(1:4)=true;                        % extract east, west, north, south
J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bcIxVE = Gt.cells.faces(I & J, 1);

% hydrostatic pressure conditions for open boundary faces
p_bc     = Gt.faces.z(bcIxVE) * rhow * g(3);
bc2D     = addBC([], bcIxVE, 'pressure', p_bc); 
bc2D.sat = repmat([1 0], numel(bcIxVE), 1);

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W2D, 'bc', bc2D);
schedule.control(2) = struct('W', W2D, 'bc', bc2D);
schedule.control(2).W.val = 0;
% Specifying length of simulation timesteps
injDuration         = 30;       % injection duration in injTimestepLength
injTimestepLength   = year;     % injection timestep length
redDuration         = 97;       % redistribution duration in injTimestepLength
redTimestepLength   = 10*year;  % redistribution timestep length

% Specifying length of simulation timesteps
schedule.step.val = [repmat(injTimestepLength, injDuration, 1); ... 
                     repmat(redTimestepLength,redDuration, 1)];
% The injection timesteps will use control 1, the redistribution timesteps 
% will use control 2.
schedule.step.control = [ones(injDuration, 1); ...
                         ones(redDuration, 1) * 2];

%% Create and simulate model
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
% setup fault block transmissibilities
% reconstructing vertically integerated transmissibilities
rock_tmp            = rock2D; 
rock_tmp.perm       = rock2D.perm .* Gt.cells.H; 
T1                  = computeTrans(Gt, rock_tmp); 
cf                  = Gt.cells.faces(:, 1); 
nf                  = Gt.faces.num; 
T1                  = 1 ./ accumarray(cf, 1 ./ T1, [nf, 1]);
setToZero           = false(Gt.faces.num, 1); 
setToZero(fblock)   = true;
T1(setToZero)       = 0;
model.operators.T   = T1(model.operators.internalConn);
model.operators.T_all(model.operators.internalConn) = model.operators.T;

% fault properties
% refer manuscript for equations and term explanations
faultpermeability       = 0.001 * milli * darcy;    % fault permeability [m2]
caprockl                = 500;                      % caprock length [m]
faultwidth              = 5;                        % fault thickness (damage zone thickness) [m] 
faultbreadth            = 200;                      % fault breath in [m]
pefault                 = 0.2*barsa;                % fault capillary entry pressure [Pa]
qleaktran = (faultwidth*faultbreadth*faultpermeability/caprockl);

%fault and caprock properties
% fault transmissibility
faulttran               = zeros(model.G.cells.num,1); 
faulttran(cblock)       = qleaktran; 
model.rock.ff           = faulttran;
% fault caprock length
faultcaprockL           = zeros(model.G.cells.num,1); 
faultcaprockL(cblock)   = caprockl;
model.rock.lc           = faultcaprockL;
% fault capillary entry pressure
faultPe                 = zeros(model.G.cells.num,1); 
faultPe(cblock)         = pefault; 
model.rock.fpe          = faultPe;
% numerical constraints - Change this parameter to prevent oscillatory
% leakage response
hlim                    = zeros(model.G.cells.num,1); 
hlim(cblock)            = 1e-6;
model.rock.hlimit       = hlim;

% simulator
ToutG = getFaceTransmissibility(Gt, rock2D);
[wellsol, states] = simulateScheduleAD(initState, model, schedule);
%% leakage post processing
time                = cumsum([schedule.step.val])/year; % time array in years
tval                = length(states);   % length of the time array (total timesteps)
tf                  = 86400*365;        % converting seconds to years
mf                  = 1e9;              % converting kg to Mt
lflux               = zeros(tval,1);    % array for fault gas leakage store [Mt/year]
influx              = zeros(tval,1);    % array for gas injection store [Mt/year]
boundfluxg          = zeros(tval,1);    % array to store gas flux through boundaries [Mt/year]
boundfluxw          = zeros(tval,1);    % array to store water flux through boundaries [Mt/year]
totalfaultleakage   = 0.0;              % cumulative fault gas leakage [Mt]    
totalwaterleakage   = 0.0;              % cumulative fault water leakage [Mt]
totalinjected       = 0.0;              % cumulative gas injected [Mt]
totalboundleakedg   = 0.0;              % cumulative gas flow through boundaries [Mt]
totalboundleakedw   = 0.0;              % cumulative water flow through boundaries [Mt]

for i = 1:tval
    lflux(i) = sum(states{i}.fltflux(cblock))*tf*co2_rho/mf;
    influx(i) = wellsol{i}.flux(2)*tf*co2_rho/mf;
    boundfluxg(i) = sum(states{i}.flux(bc2D.face,2))*tf*co2_rho/mf;
    boundfluxw(i) = sum(states{i}.flux(bc2D.face,1))*tf*rhow/mf;
    if i == 1
        totalfaultleakage   = lflux(i)*time(i);
        totalinjected       = influx(i)*time(i);
        totalboundleakedg   = boundfluxg(i)*time(i);
        totalboundleakedw   = boundfluxw(i)*time(i);
    else
        totalfaultleakage   = totalfaultleakage + lflux(i)*(time(i)-time(i-1));
        totalinjected       = totalinjected + influx(i)*(time(i)-time(i-1));
        totalboundleakedg   = totalboundleakedg + boundfluxg(i)*(time(i)-time(i-1));
        totalboundleakedw   = totalboundleakedw + boundfluxw(i)*(time(i)-time(i-1));
    end
end
% creating the new saturation term in states.
for i = 1:tval
    [h, h_max] = upscaledSat2height(states{i}.s(:,2), states{i}.sGmax, Gt, ...
                                    'pcWG', fluid.pcWG, ...
                                    'rhoW', fluid.rhoW, ...
                                    'rhoG', fluid.rhoG, ...
                                    'p', states{i}.pressure);

    states{i}.satgasfine = height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid);
    states{i}.hco2 = h;
end

%% results analysis
% grid and saturation plot
% This plot features 4 subplot. subplot 1 contains the upscaled 
% CO2 saturation at the end of year 1. Sublot 2 contains the upscaled 
% CO2 saturation at the end of year 30. Subplot 3 contains the upscaled CO2
% saturation at the end of year 1000. Subplot 4 contains the gas injection
% and fault gas leakage as a function of time.
tol = 1e-4;
figure; 
set(gcf, 'Position',  [100, 100, 1200, 800], 'Color', 'w');
subplot(2,2,1);
sat1 = states{1}.s(:,2); sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; 
clim([0.0 0.05]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off; 
title("Time = " + num2str(time(1))+" years ");
subplot(2,2,2);
sat1 = states{30}.s(:,2); sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10);
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; clim([0.0 0.05]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off; 
title("Time = " + num2str(time(30))+" years ");
subplot(2,2,3);
sat1 = states{end}.s(:,2); sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; clim([0.0 0.05]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off;
title("Time = " + num2str(time(end))+" years ");
subplot(2,2,4); 
names = ["Fault Leakage" "Injection"];
plot(time,lflux); 
hold on; 
plot(time,influx);
xlabel('Time [years]'); 
ylabel('Leakage Rate [MT/year]'); 
legend(names,'Location','best');
set(gca, 'YScale', 'log');
set(gca,'Xscale','log'); 
ylim([10^-6 10]);
grid on;
