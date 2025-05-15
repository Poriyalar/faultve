function [problem, state] = equationsWGVEbasic(model, state0, state, dt, drivingForces, varargin)
%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   opt = struct('reverseMode', false, ...
                'resOnly', false, ...
                'iteration', -1,...
                'adjointForm',false);

   opt = merge_options(opt, varargin{:});

   assert(isempty(drivingForces.src)); % unsupported

   s  = model.operators;
   f  = model.fluid;
   G  = model.G;

   % Extract the current and previous values of all variables to solve for
   if(opt.adjointForm)
      % use sGmax as primary variable
      [p, sG, sGmax, wellSol] = model.getProps(state , 'pressure', 'sg', 'sGmax', 'wellsol');
   else
      [p, sG, wellSol] = model.getProps(state , 'pressure', 'sg', 'wellsol');
   end
   [p0, sG0, sGmax0, wellSol0] = model.getProps(state0, 'pressure', 'sg','sGmax', 'wellsol');

   
   [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);


   % Initialization of independent variables

   if ~opt.resOnly
      % ADI variables needed since we are not only computing residuals
      if ~opt.reverseMode
         if(opt.adjointForm)
            [p, sG, sGmax, wellVars{:}] = initVariablesADI(p, sG, sGmax, wellVars{:});
         else
            [p, sG, wellVars{:}] = initVariablesADI(p, sG, wellVars{:});
            sGmax = max(sG, sGmax0);
         end
      else
          wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
         if(opt.adjointForm)
            [p0, sG0, sGmax0, wellVars0{:}] = initVariablesADI(p0, sG0, sGmax0, wellVars0{:});
         else
             warning('using non adjoint form for backward simulation: Hysteresis may be not properly handled')
            [p0, sG0, wellVars0{:}] = initVariablesADI(p0, sG0, wellVars0{:});
            sGmax = max(sG, sGmax0);
         end
      end
   elseif ~opt.adjointForm
      % sGmax was not yet defined.  We are probably in the process of
      % evaluating residuals after non-convergence.  Set it here to avoid
      % problems further down
      sGmax = max(sG, sGmax0);
   end
   

   sW  = 1 - sG;  % for ease of reading, we define an explicit variable
   sW0 = 1 - sG0; % also for water saturation

   % Preparing various necessary, intermediate values

   % multiplier for mobilities
   [pvMult, transMult, mobMult, pvMult0] = getMultipliers(f, p, p0);

   % relative permeability
   [krW, krG] = model.evaluateRelPerm({sW, sG}, p, 'sGmax', sGmax);
   krW = krW * mobMult;
   krG = krG * mobMult;

   % Transmissibility
   trans = s.T .* transMult;

   % Gravity gradient per face, including the gravity component resulting
   % from caprock geometry
   gdz = model.getGravityGradient();

   % CO2 phase pressure
   pW  = p;
   pW0 = p0;
   pG  = p + f.pcWG(sG, p, 'sGmax', sGmax);

   % Evaluate water and CO2 properties
   [vW, bW, mobW, rhoW, upcw] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krW, trans, gdz, 'W', 0, 0);
   [vG, bG, mobG, rhoG, upcg] = ...
       getPhaseFluxAndProps_WGVE(model, pW, pG, krG, trans, gdz, 'G', 0, 0);
   bW0 = f.bW(pW0);
   bG0 = f.bG(pW0);             % Yes, using water pressure also for gas here

   if model.outputFluxes
       state = model.storeFluxes(state, vW, [], vG);
   end
   
   % Multiply upstream b-factors by interface fluxes to obtain fluxes at
   % standard conditions
   bWvW = s.faceUpstr(upcw, bW) .* vW;
   bGvG = s.faceUpstr(upcg, bG) .* vG;


   % Setting up brine and CO2 equations

   % Water (Brine)
   eqs{1} = (s.pv / dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + s.Div(bWvW);

   % Gas (CO2)
   eqs{2} = (s.pv / dt) .* (pvMult .* bG .* sG - pvMult0 .* bG0 .* sG0) + s.Div(bGvG);

   % Setting names of variables and equations
   primaryVars = {'pressure' , 'sG'}; 
   types = {'cell' , 'cell'};
   names = {'water', 'gas'};  
   % Include influence of boundary conditions
   if state0.s(end, 2) > 0.1
       %keyboard;
   end
   
   
   %---------------------------------NewFaultModel---------------------------------

   % fault leakage model by Hariharan
   g = norm(model.gravity);                % defining gravity
   rr = model.rock;                        % rock model abbreviation
   muG = f.muG(pW);                        % gas viscosity
   bouyT = (rhoW-rhoG) .* g;               % bouyancy pressure contribution
   viscT = (1./muG);                       % viscous contribution
   gh = G.cells.H;                         % grid block height
   modelz = G.cells.z;                     % model gridblock depth
   pat = pW - rhoW .* g .* modelz;         % aqueous phase overpressure simple-ve
   % fault properties
   lc = rr.lc;                             % faulted caprock height
   pe = rr.fpe;                            % fault capillary entry pressure
   swr = f.res_water;                      % residual water saturation
   sgr = f.res_gas;                        % residual gas saturation
   hlim = rr.hlimit;                       % hgas height check for numerical stability
   plim = 1e-1;                            % pressure check for numerical stability

   pen = pe;
   knfac = 1;
   krfactor = 1;                                                        % factor to acount for pc fault effect
   
   % new method based on reconstruction - co2 height estimation
   hgas = gh .* ((sG .* (1-swr) - (sGmax .* sgr)) ./ ((1-swr) .* (1 - swr - sgr)));
   hmax = gh .* (sGmax ./(1-swr));
  
   %[hgas, hmax] = upscaledSat2height(sG,sGmax,G,'pcWG',f.pcWG,'rhoW',rhoW,'rhoG',rhoG,'p',pW);

   % fault flow flux estimation
   x1 = hgas > hlim;                                                    % leakage opens once the gas is saturated for 1micrometer with fault block
   x2 = (pat + bouyT .* hgas - pen) > plim;                             % leakage open once capillary entry pressure is breached
   krf = x1 .* x2;                                                      % dirac style switch to open leakage flux
   %krsmooth = 1.0 - exp(-1e1 .* (pat + bouyT .* hgas - pen) .^ 1.0);   % smoothening parameter
   krsmooth = 1.0 - exp(-1e-3 .* (pat + bouyT .* hgas - pen) .* x2);
   %krsmooth = 1.0;                                     % smoothening parameter
   ff = rr.ff .* krf .* krfactor .* knfac .* krsmooth;  % the static transmissibility terms for equations
   lgflux = ff .* viscT .*(pat + hgas .* bouyT + lc .* bouyT - pen);         % gas leakage flux
   
   
   
   % add equations if the fault is active
   if ~isempty(rr.ff)
       eqs{2} = eqs{2} + lgflux;
       if model.outputFluxes
       state = model.storefaultfluxes(state,lgflux);
       end
   end 

   %-----------------------------------------------------------------------
   
   [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                  {pW, pG}, {sW, sG}, {mobW, mobG}, {rhoW, rhoG}, ...
                                                  {}, {}, drivingForces);

   if(opt.adjointForm)
       eqs{3}=sGmax-max(sG,sGmax0);
   end
       


   % Setting up well equations
   dissolved = {{[],[]},{[],[]}};  % two phases, no dissolution
   %add hysteresis variable, equation and name
   if(opt.adjointForm)
       primaryVars = {primaryVars{:}, 'sGmax'}; %#ok
       types = {types{:} , 'cell' };  %#ok
       names = {names{:} , 'gasMax'}; %#ok
   end
   
   % add names of well equations
   primaryVars = {primaryVars{:}, wellVarNames{:}}; %#ok
   
   [eqs, names, types, state.wellSol] = ...
       model.insertWellEquations(eqs, names, types, wellSol0, wellSol, ...
       wellVars, wellMap, p, ...
       {mobW, mobG}, {rhoW, rhoG}, dissolved, ...
       {}, dt, opt);

   % Setting up problem
   problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
