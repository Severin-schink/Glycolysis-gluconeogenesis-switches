function r_0 = Parameters_flux_balance(phi1,phi2,efficiency)
% Compute_parameters 
% Get steady state fluxes and steady state concentrations (S0)
% Steady state fluxes are lumped into Vmax, Km based on BRENDA DB, see SI.

close all;

%% GLYCOLYTIC STEADY STATE FLUXES (r_0)and STEADY STATE Concentrations (S0)

% Biomass flux
r_BM = 0.62; % (mM/s)
r_0(12) = r_BM;

% Biomass precursor composition
beta = [0.209 0.292 1.289];

% Drain to avoid glycolytic carbon buildup
r_0(11) = 0.1*r_BM; 
% Drain to avoid TCA carbon buildup
r_0(6) = 0.02*r_BM; 

% glucose uptake
r_0(1) = r_BM + r_0(11) + 0.5*r_0(6);

% vmax acetate uptake 
r_0(9) = 0.9;

% --- Fluxes through central metabolism ---

% Net flux through the upper irreversible reaction
r_up = r_0(1) - r_0(11) - beta(1)*r_BM;

% Upper glycolysis
r_0(2) = (1+phi1)*r_up;
% Lower gluconeogenesis
r_0(3) = phi1*r_up;

%super enzyme ENO 
r_0(7) = r_up*(1/efficiency-1);    %reverse
r_0(8) = r_up/efficiency;          %forward

% After Drain through PEP for BM
r_low = r_up - 0.5*beta(2)*r_BM;

% Lower glycolysis
r_0(4) = 2*(1+phi2)*r_low; %2, due to going to 2 trioses
% Lower gluconeogenesis
r_0(5) = 2*phi2*r_low;  %2, due to going to 2 trioses, what remains after the futile cycle

end

