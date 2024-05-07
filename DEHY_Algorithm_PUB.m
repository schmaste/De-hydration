% Pseudo-transient finite difference algorithm for 1D HMC (de)hydration simulation.
% This algorithm was used to generate results for the manuscript titled:
% "(De)hydration Front Propagation into Zero-Permeability Rock"
% Code authors: Lyudmila Khakimova, Yuri Podladchikov, Stefan Schmalholz (November 2023; update May 2024)
clear variables, close all, clf,
nx              = 201;      % Number of grid points
% Simulation for dehydration (icase = 1) and hydration (icase = 2)
for icase = 1:2
    % Independent scales
    Pmax        = 1;        % Magnitude of pressure perturbation
    k_etaf0     = 1;        % Permeability divided by fluid viscosity
    t_pert      = 1;        % Duration of pressure perturbation
    % Nondimensional parameters
    rhos_rhof   = 2;        % Ratio of solid to fluid density
    Lc_Lx       = 1;        % Ratio of characteristic length to model width
    % Model parameters
    npow        = 3;        % Porosity exponent in porosity-permeability relation
    dP_Pmax     = 1/10;     % Ratio of pressure interval across the reaction to pressure perturbation
    drhosCim_r  = (-1)^icase*0.4;   % Change in solid density times mass fraction across reaction
    drhos_r     = (drhosCim_r>0)*drhosCim_r*1.1 - (drhosCim_r<0)*drhosCim_r/4; % Making sure that total density increases with pressure
    betaphiPmax = 0.01;     % Ratio of pore compressibility  to fluid pressure
    betasPmax   = 0.0;      % Ratio of solid compressibility to fluid pressure
    betafPmax   = 0.1;      % Ratio of fluid compressibility to fluid pressure
    ttime       = 10*t_pert;                    % Simulation time
    Lx          = 6*sqrt(t_pert*k_etaf0*Pmax);  % Model width
    Lc          = Lc_Lx*Lx;                     % Characteristic length
    eta_phi     = Lc^2/k_etaf0;                 % Compaction viscosity
    beta_s      = betasPmax/Pmax;               % Solid compressibility
    beta_f      = betafPmax/Pmax;               % Fluid compressibility
    dP          = sign(drhosCim_r)*dP_Pmax*Pmax;% Pressure interval across reaction
    P_ini       = 0;                            % Initial, or ambient, pressure
    Pr          = 0;                            % Reaction pressure
    drhof_r     = 0;                            % Change in fluid density across reaction
    rhosCim_r0  = 2;                            % Initial product of solid density and mass fraction
    rhof_r      = 1;                            % Fluid density
    rhos_r      = rhos_rhof*rhof_r;             % Solid density
    % Numerical parameters
    nout        = 2;                            % Steps between visualization 
    niter       = 1e7;                          % Number of maximal iterations
    tol         = 1e-3;                         % Tolerance, or error
    dP_lut      = 0.01*Pmax;                    % Pressure array-step used for parameterization
    CFL         = 1/2;                          % Courant-Friedrich-Lax criterion value
    % preprocessing
    dx          = Lx/(nx-1);                    % Grid spacing
    x           = 0:dx:Lx;                      % Coordinate vector
    if sign(drhosCim_r) > 0
        Pf_lut   =  P_ini:dP_lut:Pmax;          % Pressure vector for initial parameterization; hydration
    else
        Pf_lut   = -Pmax:dP_lut:P_ini;          % Pressure vector for initial parameterization; dehydration
    end
    % Parameterizations
    rhosCim_lut = lut(Pf_lut,beta_s,rhosCim_r0,drhosCim_r,dP,Pr);
    rhophisCim0 = min(rhosCim_lut);
    phi_lut     = 1 - rhophisCim0./rhosCim_lut;
    rhof_lut    = lut(Pf_lut,beta_f  ,rhof_r,drhof_r,dP,Pr);
    rhos_lut    = lut(Pf_lut,beta_s  ,rhos_r,drhos_r,dP,Pr);
    rhot_lut    = phi_lut.*rhof_lut + (1-phi_lut).*rhos_lut;
    phi_max     = max(phi_lut);
    beta_tot    = min(diff(rhot_lut)./diff(Pf_lut)/rhof_lut(1));
    if beta_tot < 0, disp(beta_tot);break,end % Stop if compressibilty of system is negative; not thermodynamically stable
    % Initial conditions
    Pf          = P_ini + 0*x;              % Fluid pressure
    Pf(1)       = Pmax*sign(drhosCim_r);    % Fluid pressure at left model boundary
    Pt          = P_ini + 0*x;              % Total pressure
    Vol_Vol0    = 1     + 0*x;              % Volume change
    Vs          = zeros(1,nx+1);            % Solid velocity
    time        = 0;                        % Model time
    it          = 0;                        % Number of time step
    % Time loop 
    while time < ttime
        it          = it + 1;
        % Iteration loop with numerical calculations
        for iter = 1:niter
            rhosCim = lut(Pf,beta_s,rhosCim_r0,drhosCim_r,dP,Pr);   % Solid density times mass fraction
            phi     = 1 - rhophisCim0./rhosCim./Vol_Vol0;           % Porosity
            rhof    = lut(Pf,beta_f,rhof_r,drhof_r,dP,Pr);          % Fluid density
            rhos    = lut(Pf,beta_s,rhos_r,drhos_r,dP,Pr);          % Solid density
            rhot    = phi.*rhof + (1-phi).*rhos;                    % Total density
            if iter == 1,rhot_old = rhot; rhos_old = rhos;end
            k_etaf  = k_etaf0.*avx(phi/phi_max).^npow;              % Permeability / fluid viscosity
            dt      = 1e3*CFL*dx^2/(max(k_etaf)/beta_tot);          % Adaptive time step
            if it==1;dt = 1e1*CFL*dx^2/(max(k_etaf)/beta_tot);end
            dtiter  = CFL*dx^2/(max(k_etaf)/beta_tot)*4;            % Pseudo-transient step
            qD      = -k_etaf.*(diff(Pf)/dx);                       % Fluid flux
            ResPf   = -(rhot(2:end-1) - rhot_old(2:end-1))/dt ...   % Mass conservation equation; automated up-downwind scheme
                -     diff( avx(rhof).*qD + rhot(1:end-1).*max(0,Vs(2:end-1)) + rhot(2:end).*min(0,Vs(2:end-1)) )/dx;
            Pf(2:end-1) = Pf(2:end-1) + dtiter*ResPf;               % Fluid pressure update
            if max(abs(ResPf))<tol && iter>10,break,end
        end
        divVs       = (Pf-Pt)/eta_phi;                              % Divergence of solid velocity
        Vs          = [0 cumsum(divVs*dx)];                         % Calculation of solid velocity
        Vol_Vol0    = exp(log(Vol_Vol0) + dt*divVs);                % Volume change
        time        = time + dt;                                    % Time update
        if time >= t_pert, Pf(1) = P_ini; end                       % Remove Pf-perturbation when time > t_pert
        if time > t_pert && time > 10*t_pert, break;end             % Stop simulation when time > 10*t_pert 
        if mod(it,nout) == 0 || it==1                               % Visualization
            subplot(311),plot(x,Pf,'k'),ylabel('P_f (dimensionless)'), title('Fluid pressure'), axis([0 6 -1.05 1.05])
            text(4,0.6,['t / t_c = ',num2str(time,3)],'fontsize',14)
            subplot(312),plot(x,phi,'k'),ylabel('\phi'),title('Porosity'), axis([0 6 0 0.2])
            subplot(313),plot(x,Vol_Vol0,'k'),ylabel('V / V_0'),title('Volume change'), axis([0 6 0.94 1.05])
            xlabel('x / L_c')
            if icase==1, sgtitle('DEHYDRATION'), else, sgtitle('HYDRATION'), end, drawnow
        end
    end
    pause(4)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = lut(P,beta_phi,phi_r0,dphi_r,dP,Pr) % Function for parameterization
phi = phi_r0 + (P-Pr)*beta_phi + abs(dphi_r)*(P-Pr)./(dP+P-Pr);
end
function A_av = avx(A)
A_av = 0.5*(A(1:end-1)+A(2:end));
end