%% Sloping Aquifer with leaky Fault
% The case presented is designed to test the fault leakage function 
% for a gently sloping reservoir with a large fault.
% The reservoir is discretized as a VE grid with 2500 cells,
% each measuring 20m x 20m in the x and y directions respectively, 
% resulting in a physical model with dimensions of 1000 mx 1000m. 
% The reservoir has a thickness of 10m and a slope of ~3degrees.
% For the simulation, we use the fully-implicit solver in MRST, 
% based on automatic differentiation.

% This example is described in more detail in Section 3.3 of the manuscript
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

%% Creating the rock and grid
lx = 1000; ly = 1000; lz = 10;  % reservoir lengths in [m] 
nx = 100; ny = 100; nz = 10;    % number of cells 
% the next few lines ensure that all gridblocks contain the reservoir depth
% with a slope of 3 degrees.
Reservoirdepth = 1000;      % the reservoir depth in [m] at the top of the reservoir
physdim = [lx, ly, lz]; celldim = [nx, ny, nz];
x = linspace(0, physdim(1), celldim(1)+1);
y = linspace(0, physdim(2), celldim(2)+1);
depthspace = repmat(Reservoirdepth,numel(x),numel(y));
for i = 1:numel(x)
    depthspace(i,:) = Reservoirdepth - i*(lx/nx/20);
end
G = cartGrid(celldim, physdim,'depthz', depthspace);
G = computeGeometry(G);

% The reservoir rock is created with a 100 [md] permeability and 0.3 [-] porosity
rock = makeRock(G, 100*milli*darcy, 0.3);

% top surface grid creation (VE 2D grid)
[Gt, ~, transMult] = topSurfaceGrid(G);
rock2D = averageRock(rock, Gt);

%% creating a fault block
% The gridblocks  and faces containing the fault block are identified 
% fblock contains the faces that need to be blocked due to fault core
% cblock contains the cells that hold the fault leakage terms
fblock = addSealingFacesX(Gt,'x_range',[500, 500],'y_range',[100, 900]);
cblock = find(Gt.cells.centroids(1:Gt.cells.num,1) >= (500-0.75*lx/nx) & ...
    Gt.cells.centroids(1:Gt.cells.num,1) <= 500 & ...
    Gt.cells.centroids(1:Gt.cells.num,2) >= 100 & ...
    Gt.cells.centroids(1:Gt.cells.num,2) < 900);

%% creating fluid and rock fluid parameters (reservoir)
g       = gravity;
rhow    = 1001;                                 % water density in [kg m^-3]
co2     = CO2props();                           % CO2 property functions
p_ref   = 30 *mega;                             % reference pressure in [Pa]
t_ref   = 94+273.15;                            % reference temperature in [K]
co2_rho = 389.7;                                % CO2 density in [kg m^-3]
co2_c   = 0.0;                                  % CO2 compressibility in [Pa^-1]
wat_c   = 0;                                    % water compressibility in [Pa^-1]
c_rock  = 4.35e-5 / barsa;                      % rock compressibility in [Pa^-1]
srw     = 0.27;                                 % residual water [-]
src     = 0.2;                                  % residual CO2 [-]
muw     = 3.13e-4;                              % brine viscosity in [Pa sec] 
muco2   = 3.21e-5;                              % CO2 viscosity in [Pa sec]
pe      = 0.00;                                 % capillary entry pressure in [Pa]
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^1;     % uses CO2 saturation

% create VE fluid structure
fluid   = makeVEFluid(Gt, rock, 'sharp interface'     , ...
               'co2_mu_ref'  , muco2                  , ...
               'wat_mu_ref'  , muw                    , ...
               'co2_rho_ref' , co2_rho                , ...
               'wat_rho_ref' , rhow                   , ...
               'co2_rho_pvt' , [co2_c, p_ref]         , ...
               'wat_rho_pvt' , [wat_c, p_ref]         , ...
               'residual'    , [srw, src]             , ...
               'pvMult_p_ref', p_ref                  , ...
               'pvMult_fac'  , c_rock                 , ...
               'invPc3D'     , invPc3D                , ...
               'kr3D'        , kr3D                   , ...
               'transMult'   , transMult);

%% creating well, initial state, boundary state and schedule
initState.pressure = rhow * g(3) * Gt.cells.z;
initState.s        = repmat([1, 0], Gt.cells.num, 1);
initState.sGmax    = initState.s(:,2);

inj_rate    = 1e-4;                 % injection rate in [m3/s]       
wellCell    = findEnclosingCell(Gt, [100,500]);  
W           = [];
W           = addWellVE(W, Gt, rock2D, wellCell, 'name', 'injector', ...
                'type', 'rate', ... % inject at constant rate
                'val', inj_rate, ...% volumetric injection rate
                'comp_i', [0 1]);   % inject CO2, not water)
W2D         = W;

% Find bounding indices in 2D grid
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
j(1)=true;                          % extract east, west, north, south
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

injDuration         = 3;        % Injection duration in injTimestepLength
injTimestepLength   = year;     % Injection timestep length
redDuration         = 97;       % Redistribution duration in injTimestepLength
redTimestepLength   = year;     % Redistribution timestep length

% Specifying length of simulation timesteps
schedule.step.val     = [repmat(injTimestepLength, injDuration, 1); ... 
                        repmat(redTimestepLength,redDuration, 1)];
% Specifying which control to use for each timestep.
schedule.step.control = [ones(injDuration, 1); ...
                         ones(redDuration, 1) * 2];

%% creating figure for  grid validation

figure(1); 
set(gcf, 'Position',  [100, 100, 800, 800], 'Color', 'w');
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotCellData(Gt,Gt.cells.z,'EdgeColor','None');
plotWell(Gt.parent, W,'FontSize',10);
plotFaces(Gt,fblock,'r'); 
h3 = colorbar; 
clim([950 1000]);
ylabel(h3, 'Model Depth');
xlabel('X - Distance [m]'); 
ylabel('Y - Distance [m]'); 
zlabel('Z - Distance [m]');
title("Simulation Grid");
view(3);

%% fault based simulation model
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
faultpermeability       = 0.001*milli*darcy;    % fault permeability [m2]
caprockl                = 10;                   % caprock length [m]
faultwidth              = 5;                    % fault thickness [m]
faultbreadth            = ly/ny;                % fault breath [m] set to cell size
pefault                 = 0.05*barsa;           % fault capillary entry pressure [Pa]
qleaktran               = (faultwidth*faultbreadth*faultpermeability/caprockl);
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
% numerical constraints
hlim                    = zeros(model.G.cells.num,1); 
hlim(cblock)            = 1e-6;
model.rock.hlimit       = hlim;

% simulator
ToutG = getFaceTransmissibility(Gt, rock2D);
[wellsol, states] = simulateScheduleAD(initState, model, schedule);

%% output
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
  
% storing output in array for quicker access
output(1,1) = totalinjected;
output(1,2) = totalfaultleakage;
output(1,3) = totalboundleakedg;
output(1,4) = totalboundleakedw;
output(1,5) = totalfaultleakage*100/totalinjected;

% saturation evolution
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
%% creating figure for  grid validation
% this plot features 3 subplot. subplot 1 contains the simulation grid
% subplot 2 contains the upscaled CO2 saturation at the end of injection
% subplot 3 contains the upscaled CO2 saturation at the end of simulation
tol = 1e-4;
climmax = 0.2;
figure(2); 
set(gcf, 'Position',  [100, 100, 1500, 400], 'Color', 'w');
subplot(1,3,1);
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotCellData(Gt,Gt.cells.z,'EdgeColor','None');
plotWell(Gt.parent, W,'FontSize',10);
plotFaces(Gt,fblock,'r'); 
h3 = colorbar; 
clim([950 1000]);
ylabel(h3, 'Model Depth');
xlabel('X - Distance [m]'); 
ylabel('Y - Distance [m]'); 
zlabel('Z - Distance [m]');
title("Simulation Grid"); 
view(3);
subplot(1,3,2);
sat1 = states{3}.s(:,2); 
sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; 
clim([0.0 climmax]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off; 
title("Time = " + num2str(time(3))+" years ");
subplot(1,3,3);
sat1 = states{100}.s(:,2); 
sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; 
clim([0.0 climmax]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off; 
title("Time = " + num2str(time(100))+" years ");

%% Rate plots
% this plot contains the injection and leakage rates as a function of time
names = ["Fault Gas Leakage" "Injection" "Boundary Leakage Gas" "Boundary Leakage Water"];
figure(3); 
set(gcf, 'Position',  [100, 100, 500, 500], 'Color', 'w');
plot(time,lflux); 
hold on; 
plot(time,influx); 
plot(time,boundfluxg); 
plot(time,boundfluxw);
xlabel('Time [years]'); 
ylabel('Rate [Mt/year]'); 
legend(names,'Location','best');
set(gca, 'YScale', 'log');
set(gca,'Xscale','log'); 
ylim([10^-10 10]);
grid on;

%% grid and saturation plot
% This plot features 4 subplot. subplot 1 contains the upscaled 
% CO2 saturation at the end of year 1. Sublot 2 contains the upscaled 
% CO2 saturation at the end of year 30. Subplot 3 contains the upscaled CO2
% saturation at the end of year 100. Subplot 4 contains the gas injection
% and fault gas leakage as a function of time.

tol = 1e-4;
figure(4); 
set(gcf, 'Position',  [100, 100, 1200, 800], 'Color', 'w');
subplot(2,2,1);
sat1 = states{1}.s(:,2); 
sat1(sat1<tol) = NaN;
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
sat1 = states{30}.s(:,2); 
sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; 
clim([0.0 0.05]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off;
title("Time = " + num2str(time(30))+" years ");
subplot(2,2,3);
sat1 = states{end}.s(:,2);
sat1(sat1<tol) = NaN;
plotGrid(Gt,'FaceColor','none', 'EdgeAlpha', 0.1);
plotWell(Gt.parent, W,'FontSize',10); 
plotCellData(Gt,sat1,'EdgeColor','None');
plotFaces(Gt,fblock,'r'); 
view(3); 
h2 = colorbar; 
clim([0.0 0.05]);
ylabel(h2, 'Upscaled CO2 Saturation'), axis tight off; 
title("Time = " + num2str(time(end))+" years ");
subplot(2,2,4); 
names = ["Fault Leakage" "Injection"];
plot(time,lflux); 
hold on; 
plot(time,influx);
xlabel('Time [years]'); 
ylabel('Leakage Rate [Mt/year]'); 
legend(names,'Location','best');
set(gca, 'YScale', 'log');
set(gca, 'Xscale','log'); 
ylim([10^-6 10]);
grid on;
