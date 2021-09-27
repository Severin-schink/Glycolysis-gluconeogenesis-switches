function[Modelv0,simData] = Model(Vmax,KmN,S0,Snorm,alphas,GNGconstants,GLYconstants,shifts)

%% Description of variables
% Vmax: numbering as it is in the reactions' definitions
% Vmax: (1:10) Glc_up(pgi+) pfk fbpase pyk pps PyrBMdrain superENO- superENO+ Ace_up PyrDrainOther
% KmN: (1:5)  PFK(gly^up)  FBPase(gng^up)  PYK(gly^up)   PPS(gng^low)
% S0: (1:4)   metabolite concentrations: F6P(glycolytic carbon) FBP PEP PYR (TCA carbons)
% alphas: (1:4) allosteric effectors: FBP->pyk(gly^up), PEP->fbpase(gng^up), PEP-|pfk(gly^up), FBP-|pps(gng^low)

%% Create a model named Modelv0
Modelv0 = sbiomodel('Modelv0');

% Define simulation time
configsetObj = getconfigset(Modelv0);
limitTime = 3600*shifts.time(end);
set(configsetObj, 'StopTime', limitTime);

% Increase integration accuracy, required because metabolism much faster than protein
% dynamics.
set(configsetObj.SolverOptions, 'AbsoluteTolerance',1e-15)
set(configsetObj.SolverOptions, 'RelativeTolerance',1e-12)

%% Add species
addspecies(Modelv0, 'F6P','InitialAmount',S0(1));           % Initial concentrations of mets (Glc Steady State)
addspecies(Modelv0, 'FBP','InitialAmount',S0(2));
addspecies(Modelv0, 'PEP','InitialAmount',S0(3));
addspecies(Modelv0, 'PYR','InitialAmount',S0(4));

addspecies(Modelv0, 'pfk','InitialAmount',Vmax(2));         % Initial "concentration*k_cat" of enzymes (Glc Steady State)
addspecies(Modelv0, 'fbpase','InitialAmount',Vmax(3));
addspecies(Modelv0, 'pyk','InitialAmount',Vmax(4));
addspecies(Modelv0, 'pps','InitialAmount',Vmax(5));

%% Add parameters
addparameter(Modelv0, 'F6Pnorm', Snorm(1));                    % Initial concentrations(Glc Steady State) used to normalize
addparameter(Modelv0, 'FBPnorm', Snorm(2));
addparameter(Modelv0, 'PEPnorm', Snorm(3));
addparameter(Modelv0, 'PYRnorm', Snorm(4));

addparameter(Modelv0, 'KmPFK', KmN(1));                     % Km values
addparameter(Modelv0, 'KmFBP', KmN(2));
addparameter(Modelv0, 'KmPYK', KmN(3));
addparameter(Modelv0, 'KmPPS', KmN(4));
% addparameter(Modelv0, 'KmPDH', KmN(5));

addparameter(Modelv0, 'VmaxPDH', Vmax(6));                  % Vmax values
addparameter(Modelv0, 'VmaxR', Vmax(7));
addparameter(Modelv0, 'VmaxF', Vmax(8));
addparameter(Modelv0, 'VmaxAce', Vmax(9));

addparameter(Modelv0, 'VmaxPyrDrain', Vmax(6));            % Drain of Pyr to avoid accumulation of Pyr
addparameter(Modelv0, 'VmaxF6PDrain', Vmax(11));            % Drain of F6P to avoid accumulation of F6P

addparameter(Modelv0, 'VmaxGlc', Vmax(1));                  % Glucose uptake

addparameter(Modelv0, 'pfkPEPalpha', alphas(1));            % allosteric effector PEP -| pfk
addparameter(Modelv0, 'fbpasePEPalpha', alphas(2));         % allosteric effector PEP -> fbpase
addparameter(Modelv0, 'pykFBPalpha', alphas(3));            % allosteric effector FBP -> pyk
addparameter(Modelv0, 'ppsFBPalpha', alphas(4));            % allosteric effector FBP -| pps

addparameter(Modelv0, 'kBM', Vmax(12));

% add constants of Gluconeogenic(GNG) enzymes
addparameter(Modelv0, 'fbpaseConstant', GNGconstants(1));   % FBP0 constant for fbpase, calculated from Steady State
addparameter(Modelv0, 'ppsConstant', GNGconstants(2));      % FBP0 constant for pck, calculated from Steady State
addparameter(Modelv0, 'phiBeta_fbpase', GNGconstants(3));   % phiBeta constant(fraction of BM to this enzyme for fbpase)
addparameter(Modelv0, 'phiBeta_pps', GNGconstants(4));      % phiBeta constant(fraction of BM to this enzyme for pps)

% add constants of Glycolytic(GLY) enzymes -- fix comments below
addparameter(Modelv0, 'pfkConstant', GLYconstants(1));   % FBP0 constant for pfk, calculated from Steady State
addparameter(Modelv0, 'pykConstant', GLYconstants(2));      % FBP0 constant for pyk, calculated from Steady State
addparameter(Modelv0, 'phiBeta_pfk', GLYconstants(3));   % phiBeta constant(fraction of BM to this enzyme for pfk)
addparameter(Modelv0, 'phiBeta_pyk', GLYconstants(4));      % phiBeta constant(fraction of BM to this enzyme for pyk)

addparameter(Modelv0, 'SignalGlc', shifts.glc(1),'ConstantValue', 0);   % params-signal for implementing the carbon switch
addparameter(Modelv0, 'SignalAce', shifts.ace(1),'ConstantValue', 0);   % params-signal for implementing the carbon switch

%% Add the basic reactions to the model

% Glc uptake and drain of F6P
reactionObj = addreaction(Modelv0, 'null -> F6P', 'ReactionRate', 'VmaxGlc*SignalGlc');
reactionObj = addreaction(Modelv0, 'F6P -> null', 'ReactionRate', 'VmaxF6PDrain*(F6P/F6Pnorm)');

% PFK-FBPase junction
reactionObj = addreaction(Modelv0, 'F6P -> FBP', 'ReactionRate', 'pfk*((PEP/PEPnorm)^(-pfkPEPalpha))*(F6P/F6Pnorm)/(KmPFK + F6P/F6Pnorm)');
reactionObj = addreaction(Modelv0, 'FBP -> F6P', 'ReactionRate', 'fbpase*((PEP/PEPnorm)^(fbpasePEPalpha))*(FBP/FBPnorm)/(KmFBP + FBP/FBPnorm)');

% Reversible superENO enzyme
reactionObj = addreaction(Modelv0, 'FBP -> 2 PEP', 'ReactionRate', 'VmaxF*(FBP/FBPnorm)');
reactionObj = addreaction(Modelv0, '2 PEP -> FBP', 'ReactionRate', 'VmaxR*(PEP/PEPnorm)*(PEP/PEPnorm)');

% PYK-PPS junction
reactionObj = addreaction(Modelv0, 'PEP -> PYR', 'ReactionRate', 'pyk*((FBP/FBPnorm)^pykFBPalpha)*(PEP/PEPnorm)/(KmPYK + (PEP/PEPnorm))');
reactionObj = addreaction(Modelv0, 'PYR -> PEP', 'ReactionRate', 'pps*((FBP/FBPnorm)^(-ppsFBPalpha)*PYR/PYRnorm)/(KmPPS + (PYR/PYRnorm))');

% Acetate uptake and PYR drain (pdh)
reactionObj = addreaction(Modelv0, 'null -> PYR', 'ReactionRate', 'VmaxAce*SignalAce');
reactionObj = addreaction(Modelv0, 'PYR -> null', 'ReactionRate', 'VmaxPyrDrain*(PYR/PYRnorm)');

%% Add Biomass proxy reaction

% Since the flux from PEP and PYR is double than F6P (trioses vs. hexose) the stoichiometry is double (in order to have
% steady state)

reactionObj = addreaction(Modelv0, '0.209 F6P + 0.292 PEP + 1.289 PYR -> BM', 'ReactionRate', 'kBM*(F6P/F6Pnorm)*(PEP/PEPnorm)*(PYR/PYRnorm)');

%% Protein Expression reactions ("transcriptional" regulation)

%---- Gluconeogenic Enzymes ----%
reactionObj = addreaction(Modelv0, 'null <-> fbpase', 'ReactionRate', '0.9/3600/0.62*kBM*(F6P/F6Pnorm)*(PEP/PEPnorm)*(PYR/PYRnorm)*(phiBeta_fbpase*(1+fbpaseConstant-fbpaseConstant*(FBP/FBPnorm)) - fbpase)');
reactionObj = addreaction(Modelv0, 'null <-> pps', 'ReactionRate', '0.9/3600/0.62*(F6P/F6Pnorm)*(PEP/PEPnorm)*(PYR/PYRnorm)*(phiBeta_pps*(1+ppsConstant-ppsConstant*(FBP/FBPnorm)) - pps)');


%---- Glycolytic Enzymes ----%
reactionObj = addreaction(Modelv0, 'null <-> pfk', 'ReactionRate', '0.9/3600/0.62*(F6P/F6Pnorm)*(PEP/PEPnorm)*(PYR/PYRnorm)*(phiBeta_pfk*(pfkConstant+(1-pfkConstant)*(FBP/FBPnorm)) - pfk)');
reactionObj = addreaction(Modelv0, 'null <-> pyk', 'ReactionRate', '0.9/3600/0.62*(F6P/F6Pnorm)*(PEP/PEPnorm)*(PYR/PYRnorm)*(phiBeta_pyk*(pykConstant+(1-pykConstant)*(FBP/FBPnorm)) - pyk)');
%% Events implement carbon shift

for i = 1:length(shifts.time)-1
    addparameter(Modelv0, ['glc_' num2str(i)], shifts.glc(i)); 
    addparameter(Modelv0, ['ace_' num2str(i)], shifts.ace(i)); 
    addparameter(Modelv0, ['time_shift_' num2str(i)], shifts.time(i)); 
    addevent(Modelv0, ['time>= time_shift_' num2str(i) '*3600'], {['SignalGlc = ' ['glc_' num2str(i)]],['SignalAce = ' 'ace_' num2str(i)]});
end



% Solve model, numerical for shifts, or steadystate 
if length(shifts.time)>1
    simData = sbiosimulate(Modelv0);
    
    % Resample data to minutes
    simData = resample(simData, 0:60:limitTime);
else 

    [S,simData] = sbiosteadystate(Modelv0,'Method','simulation','AbsTol',1e-15,'RelTol',1e-12);
    
end

end

