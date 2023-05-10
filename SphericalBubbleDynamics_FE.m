clear; clc;

%%
% Parameters to be specified
%
% Undeformed geometry
R0 = 22.5e-6; % Undeformed (equilibrium) bubble radius (m)
Rout = 201*R0; % Undeformed outer radius of the simulation domain (m)

% Far-field thermodynamic conditions
P_inf = 101325; % Atmospheric pressure (Pa)
T_inf = 298.15; % Far field temperature (K)

% Type of cavitation driving force
%     LIC = Laser induced cavitation
%     UIC = Ultrasound induced cavitation
cav_type = 'LIC';

% Driving conditions
if strcmp(cav_type,'LIC') == 1
    Rmax = 4*R0; % Maximum (initial) bubble radius for LIC (m)
    %
    PA = 0; % Not used for LIC - set to be zero
    omega = 2*pi*(1e6); % Not used for LIC
    delta = pi/omega; % Not used for LIC
    n = 3.7; % Not used for LIC
elseif strcmp(cav_type,'UIC') == 1
    PA = -50*P_inf; % Amplitude of the ultrasound pulse (Pa)
    omega = 2*pi*(1e6); % Frequency of the ultrasound pulse (1/s)
    delta = pi/omega; % Time shift for the ultrasound pulse (s)
    n = 3.7; % Exponent that shapes the ultrasound pulse
    %
    Rmax = R0; % Not used for UIC
else
    disp('Incorrect cavitation type');
    return;
end

% Simulation time and output frames
tspan = 60e-6; % Total time span (s)
nframes = 400; % Number of frames in time to print progress

% Type of constitutive model structure
%     KV = Kelvin-Voigt
%     SNS = Standard nonlinear solid
mat_type = 'KV';

% Material parameters for the surrounding material
if strcmp(mat_type,'KV') == 1
    %
    %                            |
    %                      -------------
    %                      |           |
    %                      |           |
    %    Quadratic-I1      \           |
    %       spring:        /         | - | Newtonian dashpot: mu
    %                      \         |___|
    %   G,alpha,C_long     /           |
    %                      |           |
    %                      |           |
    %                      -------------
    %                            |
    %
    G = 10e3; % Ground-state shear modulus (Pa)
    alpha = 0.0; % Strain-stiffening parameter
    C_long = 1430; % Longitudinal wave speed (m/s)
    rho = 1060; % Density (kg/m^3)
    mu = 0.03; % Viscosity (Pa-s)
    %
    G1 = 10e6; % Not used for KV
    gdot0 = G/mu; % Not used for KV
    mp = 1; % Not used for KV
elseif strcmp(mat_type,'SNS') == 1
    %
    %                            |
    %                      -------------
    %                      |           |
    %                      |           \   Hencky spring:
    %    Quadratic-I1      \           /       G1
    %       spring:        /           \
    %                      \           |
    %   G,alpha,C_long     /         | - | non-Newtonian dashpot:
    %                      |         |___|     gdot0,mp
    %                      |           |
    %                      -------------
    %                            |
    %
    G = 10e3; % Ground-state shear modulus (Pa)
    alpha = 0.0; % Strain-stiffening parameter
    C_long = 1430; % Longitudinal wave speed (m/s)
    rho = 1060; % Density (kg/m^3)
    G1 = 10e6; % Non-equilibrium shear modulus (Pa)
    gdot0 = 3.3333e+05; % Reference strain rate (1/s)
    mp = 1.00; % Strain-rate sensitivity exponent
    %
    mu = G/gdot0; % Not used for SNS
else
    disp('Incorrect constitutive model type');
    return;
end

% Surface tension
gamma = 0.056; % Surface tension (N/m)

% Parameters for the bubble contents
D0 = 24.2e-6; % Binary diffusion coeff (m^2/s)
kappa = 1.4; % Specific heats ratio
Ru = 8.3144598; % Universal gas constant (J/mol-K)
Rv = Ru/(18.01528e-3); % Gas constant for vapor (Ru/molecular weight) (J/kg-K)
Ra = Ru/(28.966e-3); % Gas constant for air (Ru/molecular weight) (J/kg-K)
A = 5.28e-5; % Thermal conductivity parameter (W/m-K^2)
B = 1.17e-2; % Thermal conductivity parameter (W/m-K)
P_ref = 1.17e11; % Reference pressure (Pa)
T_ref = 5200; % Reference temperature (K)

% Numerical parameters
nElem = 1000; % Number of finite elements in the surroundings
rel_dt = 1/3; % Time increment as a fraction of the stable time increment
NT = 500; % Grid points inside the bubble, resolution should be >=500

%%
% Intermediate calculated variables
%
% General characteristic scales
if strcmp(cav_type,'LIC') == 1
    Rc = Rmax; % Characteristic length scale (m)
    Uc = sqrt(P_inf/rho); % Characteristic velocity (m/s)
    tc = Rmax/Uc; % Characteristic time scale (s)
elseif strcmp(cav_type,'UIC') == 1
    Rc = R0; % Characteristic length scale (m)
    Uc = sqrt(P_inf/rho); % Characteristic velocity (m/s)
    tc = R0/Uc; % Characteristic time scale (s)
end

% Bulk modulus
Bmod = rho*C_long^2 - (4/3)*G;

% Parameters for the bubble contents
Pv = P_ref*exp(-T_ref./T_inf); % Vapor pressure evaluated at the far field temperature (Pa)
K_inf = A*T_inf+B; % Thermal conductivity evaluated at the far field temperature (W/m-K)

%%
% Non-dimensional variables
%
C_star = C_long/Uc; % Dimensionless wave speed
We = P_inf*Rc/(2*gamma); % Weber number
Ca = P_inf/G; % Cauchy number
Reynolds = P_inf*Rc/(mu*Uc); % Reynolds number
fom = D0/(Uc*Rc); % Mass Fourier number
chi = T_inf*K_inf/(P_inf*Rc*Uc); % Lockhart–Martinelli number
A_star = A*T_inf/K_inf; % Dimensionless A parameter
B_star = B/K_inf; % Dimensionless B parameter (Note that A_star+B_star=1.)
Pv_star = Pv/P_inf; % Dimensionless vapor saturation pressure at the far field temperature
Bmod_star = Bmod/P_inf; % Dimensionless bulk modulus
G1_star = G1/P_inf; % Dimensionless non-equilibrium shear modulus
gdot0_star = gdot0*tc; % Dimensionless reference strain rate
%
tspan_star = tspan/tc; % Dimensionless time span
%
% Non-dimensional variable only used for LIC
Req = R0/Rmax; % Dimensionless equilibrium bubble radius
%
% Non-dimensional variables only used for UIC
PA_star = PA/P_inf; % Dimensionless amplitude of the ultrasound pulse
omega_star = omega*tc; % Dimensionless frequency of the ultrasound pulse
delta_star = delta/tc; % Dimensionless time shift for the ultrasound pulse

% Non-dimensional properties vector for the FE calculation in the surroundings
if strcmp(mat_type,'KV') == 1
    props = struct('G',1/Ca,'bulkmod',Bmod_star,'alpha',alpha,'shearvisco',1/Reynolds,'bulkvisco',0.003*Uc/(P_inf*Rc));
elseif strcmp(mat_type,'SNS') == 1
    props = struct('G',1/Ca,'bulkmod',Bmod_star,'alpha',alpha,'G1',G1_star,'gdot0',gdot0_star,'mp',mp,'bulkvisco',0.0);
end

% Non-dimensional parameters vector for the bubble contents
params = [NT Rv Ra kappa fom chi A_star B_star Pv_star];

%%
% Finite element spatial discretization in the surroundings
%
nNodese = 2; % Number of nodes in each element
nNodes = (nNodese-1)*nElem + 1; % Total number of nodes
%
delr = (Rout - R0)/R0; % Dimensionless domain size of the surroundings (normalized by R0)
%
% Dimensionless undeformed nodal coordinates (normalized by the relevant characteristic length scale)
%
coords = zeros(1,nNodes);
for ii = 1:nNodes
    coords(ii) = R0*(1 + delr*(ii-1)/(nNodes-1))/Rc;
end
%
% Element connectivity (specifies node numbers on each element)
%
connect = zeros(nNodese,nElem);
for ii=1:nElem
    connect(1,ii) = ii;
    connect(2,ii) = ii+1;
end

%% 
% Initial conditions for displacement, velocity, Maxwell-element internal variables, temperature, pressure, and concentration
%
% Initial condition for dimensionless displacement field
if strcmp(cav_type,'LIC') 
    u0 = (coords.^3 + 1^3 - Req^3).^(1/3)-coords; % Normalized by Rmax
elseif strcmp(cav_type,'UIC') 
    u0 = zeros(1,nNodes);
end

% Initial condition for velocity field
v0 = zeros(nNodes,1);

% Initialization of the internal variable fields for the Maxwell element (only used for mat_type='SNS')
lambdav = ones(nElem,1); % Viscous hoop stretch. Placeholder. Initialized for a fully relaxed Maxwell element when the element subroutine is first called.
Mandel = zeros(nElem,1); % Hoop component of the elastic Mandel stress. Initialized to be zero for a fully relaxed Maxwell element.

% Initial dimensionless temperature field
Theta0 = zeros(1,NT); 

% Initial bubble pressure
if strcmp(cav_type,'LIC') == 1
    P0 = Pv + (P_inf + 2*gamma/R0 - Pv)*((R0/Rmax)^3); % Initial bubble pressure for LIC
    P0_star = P0/P_inf; % Dimensionless initial bubble pressure for LIC
elseif strcmp(cav_type,'UIC') == 1
    P0 = P_inf + 2*gamma/R0; % Initial bubble pressure for UIC
    P0_star = P0/P_inf; % Dimensionless initial bubble pressure for UIC
end

% Initial vapor mass fraction field
k0 = ((1+(Rv/Ra)*(P0_star/Pv_star-1))^(-1))*ones(1,NT);

% Place the initial conditions for the bubble contents in the state vector
X0 = [P0_star Theta0 k0];
X0 = reshape(X0,length(X0),1);

% Initial traction loads for the finite element calculation
patm = 1; % Dimensionless atmospheric pressure
pin = P0_star; % Dimensionless initial pressure at the bubble wall
pout = 0; % Dimensionless initial far-field gauge pressure

%%
% Time incrementation
%
dt_star_old = rel_dt*min(diff(coords+u0))/C_star; % Initial dimensionless time increment. Taken to be rel_dt*(stable time increment).
% Make sure that the initial time step is sufficiently small for UIC
if strcmp(cav_type,'UIC') == 1
    dt_star_old = min(dt_star_old,delta_star/100);
end

% File for writing radius versus time results
if strcmp(mat_type,'KV') == 1
    fid1 = fopen(strcat(sprintf('%s_G_%1.0fkPa_alpha_%1.0f_mu_%1.0f.dat',mat_type,G/1000,alpha*100,mu*100)),'w');
elseif strcmp(mat_type,'SNS') == 1
    fid1 = fopen(strcat(sprintf('%s_G_%1.0fkPa_alpha_%1.0f_gdot0_%1.0f_mp_%1.0f.dat',mat_type,G/1000,alpha*100,gdot0,mp*100)),'w');
end
fprintf(fid1,'%.8g,%.8g,%.8g,%.8g\r\n',0,(coords(1) + u0(1)),v0(1),pin);

%% 
% Calculate the initial acceleration and midincrement velocity
%
% Allocate the global residual vector and mass matrix
%
M = zeros(nNodes,nNodes);
R = zeros(nNodes,1);
%
% Loop over the elements
%
for ii=1:nElem
    %
    % Extract the coords, displacements, and velocities of each node in the current element
    %
    coordse = zeros(1,nNodese);
    ue = zeros(1,nNodese);
    ve = zeros(1,nNodese);
    for jj = 1:nNodese
        coordse(jj) = coords(connect(jj,ii));
        ue(jj) = u0(connect(jj,ii));
        ve(jj) = v0(connect(jj,ii));
    end
    %
    % Calculate the mass matrix and the residual vector for the current element
    %
    if strcmp(mat_type,'KV') == 1
        [Me,Re] = Element_KV(nNodese,coordse,ue,ve,props);
    elseif strcmp(mat_type,'SNS') == 1
        lambdav_old = lambdav(ii);
        Me_old = Mandel(ii);
        [Me,Re,lambdav_new,Me_new] = Element_SNS(nNodese,coordse,ue,ve,props,0,lambdav_old,Me_old);
        lambdav(ii) = lambdav_new;
        Mandel(ii) = Me_new;
    end
    %
    % Add the stiffness and residual from the current element into global matrices
    % Form the global stiffness and residual matrices
    %
    for kk = 1:nNodese
        rw = connect(kk,ii);
        R(rw) = R(rw) + Re(kk);
        for ll = 1:nNodese
            cl = connect(ll,ii);
            M(rw,cl) = M(rw,cl) + Me(kk,ll);
        end
    end
    %
end
%
% Lumped mass using the row-sum method
%
Mlump = zeros(nNodes,1);
for ii = 1:nNodes
    Mlump(ii) = sum(M(ii,:));
end
%
% Add the contributions to the residual due to the inner and outer pressures
%
R(1) = R(1) + (pin - 1/(We*(coords(1) + u0(1))) - patm - pout)*(4*pi*(coords(1) + u0(1))^2);
%
% Add a contribution to the residual to damp outgoing longitudinal waves
%
R(nNodes) = R(nNodes) - C_star*v0(nNodes)*(4*pi*(coords(nNodes)+u0(nNodes))^2);
%
% Initialize the displacement, acceleration, and midincrement velocity fields
%
ut = u0'; % Initial displacement field (t=0)
at = R./Mlump; % Initial acceleration field (t=0)
vt = v0 - (dt_star_old/2)*at; % Initial midincrement velocity field (t=-dt_star/2)

%%
% Begin loop over timesteps
%
% Allocate the acceleration, velocity, and displacement vectors
%
atau = zeros(nNodes,1);
vtau = zeros(nNodes,1);
utau = zeros(nNodes,1);
%
time = 0; % Initial time
inc = 0; % Increment counter
frame = 1; % Frame counter for displaying progress
%
tic; % Start timer
%
while (time < tspan_star)
    %
    % Update time
    %
    dt_star = rel_dt*min(diff(coords+ut'))/C_star; % Dimensionless time increment. Taken to be rel_dt*(stable time increment).
    % Make sure that the time step is sufficiently small to resolve the UIC pressure pulse
    if (strcmp(cav_type,'UIC') == 1) && (time < 2*delta_star)
        dt_star = min(dt_star,delta_star/100);
    end
    time  = time + dt_star; % Dimensionless time at the end of the increment
    %
    % For UIC, calculate the time dependent far field pressure at the end of the increment
    %
    if (strcmp(cav_type,'UIC'))
        %
        if (abs(time-delta_star)>(pi/omega_star))
            pout = 0;
        elseif (abs(time-delta_star)<=(pi/omega_star))
            pout = PA_star*(0.5*(1+cos(omega_star*(time-delta_star))))^n;
        end        
        %
    end    
    %
    % Save the bubble radius at the beginning of the time increment
    %
    Rstart = coords(1) + ut(1);
    %
    % Calculate the midincrement velocity and the displacement at the end of the increment
    %
    vtau = vt + (dt_star+dt_star_old)*at/2; % At t=time-dt_star/2
    utau = ut + dt_star*vtau; % At t=time
    %
    % Save the bubble radius at the end of the time increment
    %
    Rend = coords(1) + utau(1);
    %
    %**********************************************************
    % Bubble physics
    %
    % Integrate the governing equations inside the bubble. A smaller time 
    %   increment than dt_star is used to integrate the governing equations 
    %   inside the bubble, so the bubble radius is linearly interpolated 
    %   over the course of the time increment.
    %
    % Estimate the stable time increment for the mass diffusion equation
    dt_star_mass = (min(Rstart,Rend)/NT)^2/fom;
    % Estimate the stable time increment for the temperature equation
    Theta = X0(2:(NT+1));
    T = (A_star - 1 + sqrt(1+2.*A_star.*Theta))./A_star;
    K_star = A_star.*T+B_star;
    dt_star_temp = ((min(Rstart,Rend)/NT)^2/chi)*(X0(1)*kappa/(kappa-1))/max(K_star.*T);
    % Use the smallest of dt_star, rel_dt*dt_star_mass, and rel_dt*dt_star_temp
    numSteps = ceil(dt_star/(min([rel_dt*dt_star_mass,rel_dt*dt_star_temp,dt_star])));
    dt_star_bubble = dt_star/numSteps;
    % Initialize the bubble state vector and time
    X = X0;
    t_curr = time-dt_star;
    % Step forward in time using RK4
    for ii = 1:numSteps
        k1 = bubble(t_curr,X,params,time-dt_star,time,Rstart,Rend);
        k2 = bubble(t_curr+0.5*dt_star_bubble,X+0.5*k1*dt_star_bubble,params,time-dt_star,time,Rstart,Rend);
        k3 = bubble(t_curr+0.5*dt_star_bubble,X+0.5*k2*dt_star_bubble,params,time-dt_star,time,Rstart,Rend);
        k4 = bubble(t_curr+dt_star_bubble,X+k3*dt_star_bubble,params,time-dt_star,time,Rstart,Rend);
        X = X + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dt_star_bubble;
        t_curr = t_curr+dt_star_bubble;
    end
    % Update the bubble pressure and the bubble state vector
    pin = X(1);
    X0 = X;
    %
    % Check for errors
    %
    if ~isreal(pin)
        disp('Error: Bubble pressure is imaginary');
        return;
    end
    %**********************************************************
    %
    % Allocate the global residual vector and mass matrix
    %
    M = zeros(nNodes,nNodes);
    R = zeros(nNodes,1);
    %
    % Loop over the elements
    %
    for ii=1:nElem
        %
        % Extract the coords, displacements, and velocities of each node in the current element
        %
        coordse = zeros(1,nNodese);
        ue = zeros(1,nNodese);
        ve = zeros(1,nNodese);
        for jj = 1:nNodese
            coordse(jj) = coords(connect(jj,ii));
            ue(jj) = utau(connect(jj,ii));
            ve(jj) = vtau(connect(jj,ii));
        end
        %
        % Calculate the mass matrix and the residual vector for the current element
        %
        if strcmp(mat_type,'KV') == 1
            [Me,Re] = Element_KV(nNodese,coordse,ue,ve,props);
        elseif strcmp(mat_type,'SNS') == 1
            lambdav_old = lambdav(ii);
            Me_old = Mandel(ii);
            [Me,Re,lambdav_new,Me_new] = Element_SNS(nNodese,coordse,ue,ve,props,dt_star,lambdav_old,Me_old);
            lambdav(ii) = lambdav_new;
            Mandel(ii) = Me_new;
        end
        %
        % Add the stiffness and residual from the current element into global matrices
        %
        for kk = 1:nNodese
            rw = connect(kk,ii);
            R(rw) = R(rw) + Re(kk);
            for ll = 1:nNodese
                cl = connect(ll,ii);
                M(rw,cl) = M(rw,cl) + Me(kk,ll);
            end
        end
        %
    end
    %
    % Check for errors
    %
    if ~isreal(Re)
        disp('Error: Solution is imaginary');
        return
    end
    %
    if isnan(Re)
        disp('Error: Solution is NaN');
        return
    end
    %
    % Lumped mass using the row-sum method
    %
    Mlump = zeros(nNodes,1);
    for ii = 1:nNodes
        Mlump(ii) = sum(M(ii,:));
    end
    %
    % Add the contributions to the residual due to the inner and outer pressures
    %
    R(1) = R(1) + (pin - 1/(We*Rend) - patm - pout)*(4*pi*(coords(1) + utau(1))^2);
    %
    % Add a contribution to the residual to damp outgoing longitudinal waves
    %
    R(nNodes) = R(nNodes) - C_star*vtau(nNodes)*(4*pi*(coords(nNodes)+utau(nNodes))^2);
    %
    % Update the acceleration
    %
    atau = R./Mlump;
    %
    % Save output
    %
    Uend = vtau(1)+(1/2)*dt_star*atau(1); % Velocity at the end of the increment for output
    fprintf(fid1,'%.8g,%.8g,%.8g,%.8g\r\n',time,Rend,Uend,pin);
    %
    % Prepare for the next time increment
    %
    dt_star_old = dt_star;
    ut = utau;
    vt = vtau;
    at = atau;  
    inc = inc + 1;
    %
    % Display progress
    %
    if (time >= (frame*tspan_star/nframes))
        fprintf('Output frame %i of %i\n',frame,nframes);
        fprintf('Increment = %i, Dimensionless time = %.4g, Dimensionless radius = %.4g\n',inc,time,Rend);
        toc;
        frame = frame + 1;
    end
    %
end

fclose(fid1);

%% 
% Plot the Radius vs time curve 
%
if strcmp(mat_type,'KV') == 1
    Data1 = load(strcat(sprintf('%s_G_%1.0fkPa_alpha_%1.0f_mu_%1.0f.dat',mat_type,G/1000,alpha*100,mu*100)));
elseif strcmp(mat_type,'SNS') == 1
    Data1 = load(strcat(sprintf('%s_G_%1.0fkPa_alpha_%1.0f_gdot0_%1.0f_mp_%1.0f.dat',mat_type,G/1000,alpha*100,gdot0,mp*10)));
end
%
figure(1)
plot(Data1(:,1),Data1(:,2),'k-','linewidth',2);
hold on;
set(gca,'Fontsize',12);
xlim([0 tspan_star]);
%
if strcmp(cav_type,'LIC') == 1
    xlabel('Normalized Time, $t^* = (t/R_{\rm max})\sqrt{p_\infty/\rho_0}$','Fontsize',16,'Interpreter','latex');
    ylabel('Normalized Radius, $R^*=R/R_{\rm max}$','Fontsize',16,'Interpreter','latex');
    ylim([0 1]);
elseif strcmp(cav_type,'UIC') == 1
    xlabel('Normalized Time, $t^* = (t/R_0)\sqrt{p_\infty/\rho_0}$','Fontsize',16,'Interpreter','latex');
    ylabel('Normalized Radius, $R^*=R/R_0$','Interpreter','latex','Fontsize',16);
end

%%
% Element subroutine for Kelvin-Voigt-type material models
function [Me,Re] = Element_KV(nNodese,coordse,ue,ve,props)

% Material parameters
%
shearmod = props.G;
bulkmod = props.bulkmod;
alpha = props.alpha;
shearvisco = props.shearvisco;
bulkvisco = props.bulkvisco;

% Integration points and weights
%
nInt = 1;
w = 2.;
xi = 0.;

% Calculate the current nodal coordinates
%
coordsce = coordse + ue;

% Initialize the mass matrix and the residual vector for the current element
%
Me = zeros(nNodese,nNodese);
Re = zeros(nNodese,1);

% Loop over the integration points and assemble the element stiffness matrix
%
for jj = 1:nInt
    %
    % Compute the shape function and its derivative at the current
    % integration point
    %
    sh = zeros(1,nNodese);
    dshdxi = zeros(1,nNodese);
    sh(1) = 0.5*(1.-xi(jj));
    sh(2) = 0.5*(1.+xi(jj));
    dshdxi(1) = -0.5;
    dshdxi(2) = 0.5;
    %
    % Compute dx/dxi, J, and dN/dx
    %
    dxdxi = 0.;
    dxcdxi = 0.;
    for kk = 1:nNodese
        dxdxi = dxdxi + dshdxi(kk)*coordse(kk);
        dxcdxi = dxcdxi + dshdxi(kk)*coordsce(kk);
    end
    %
    J = abs(dxdxi);
    Jc = abs(dxcdxi);
    %
    dshdx = zeros(1,nNodese);
    dshdxc = zeros(1,nNodese);
    for kk = 1:nNodese
        dshdx(kk) = dshdxi(kk)/dxdxi;
        dshdxc(kk) = dshdxi(kk)/dxcdxi;
    end
    %
    % Calculate the radial coordinate
    %
    r = 0.;
    rc = 0;
    for kk = 1:nNodese
        r = r + sh(kk)*coordse(kk);
        rc = rc + sh(kk)*coordsce(kk);
    end
    %
    % Calculate the stretches
    %
    lambdar = 1.;
    for kk = 1:nNodese
        lambdar = lambdar + dshdx(kk)*ue(kk);
    end
    %
    lambdahoop = rc/r;
    detF = lambdar*lambdahoop^2;
    %
    % Compute the spatial velocity gradient
    %
    epsrdot = 0.;
    epshoopdot = 0.;
    for kk = 1:nNodese
        epsrdot = epsrdot + dshdxc(kk)*ve(kk);
        epshoopdot = epshoopdot + sh(kk)*ve(kk)/rc;
    end
    %
    % Compute the Cauchy stress
    %
    %---------- quadratic KV model (includes bulk viscosity contribution if needed) -----------------%
    quad_fac = 1 + alpha*((lambdar^2 + 2*lambdahoop^2)/detF^(2/3) - 3 );
    sigmar = (shearmod/detF^(5/3))*quad_fac*(lambdar^2 - (lambdar^2 + 2*lambdahoop^2)/3) + bulkmod*(detF-1) ...
        + 2*shearvisco*epsrdot + (bulkvisco-(2/3)*shearvisco)*(epsrdot + 2*epshoopdot);
    sigmahoop = (shearmod/detF^(5/3))*quad_fac*(lambdahoop^2 - (lambdar^2 + 2*lambdahoop^2)/3) + bulkmod*(detF-1) ...
        + 2*shearvisco*epshoopdot + (bulkvisco-(2/3)*shearvisco)*(epsrdot + 2*epshoopdot);
    %
    % Add contribution to element residual and mass matrix from
    % current integration point
    %
    for kk = 1:nNodese
        %
        %---------- Element residual vector -----------------%
        Re(kk) = Re(kk) - sigmar*w(jj)*Jc*dshdxc(kk)*(4*pi*rc^2) - 2*sigmahoop*w(jj)*Jc*sh(kk)/rc*(4*pi*rc^2);
        %
        %------------ Element mass matrix -------------------%
        for ll = 1:nNodese
            Me(kk,ll) = Me(kk,ll) + w(jj)*J*sh(kk)*sh(ll)*(4*pi*r^2);
        end
        
    end
    %
end

end

%%
% Element subroutine for standard-nonlinear-solid-type material models
function [Me,Re,lambdav_new,Me_new] = Element_SNS(nNodese,coordse,ue,ve,props,dt,lambdav_old,Me_old)

% Material parameters
%
shearmod = props.G;
bulkmod = props.bulkmod;
alpha = props.alpha;
G1 = props.G1;
gdot0 = props.gdot0;
mp = props.mp;
bulkvisco = props.bulkvisco;

% Integration points and weights
%
nInt = 1;
w = 2.;
xi = 0.;

% Calculate the current nodal coordinates
%
coordsce = coordse + ue;

% Initialize the mass matrix and the residual vector for the current element
%
Me = zeros(nNodese,nNodese);
Re = zeros(nNodese,1);

% Loop over the integration points and assemble the element stiffness matrix
%
for jj = 1:nInt
    %
    % Compute the shape function and its derivative at the current
    % integration point
    %
    sh = zeros(1,nNodese);
    dshdxi = zeros(1,nNodese);
    sh(1) = 0.5*(1.-xi(jj));
    sh(2) = 0.5*(1.+xi(jj));
    dshdxi(1) = -0.5;
    dshdxi(2) = 0.5;
    %
    % Compute dx/dxi, J, and dN/dx
    %
    dxdxi = 0.;
    dxcdxi = 0.;
    for kk = 1:nNodese
        dxdxi = dxdxi + dshdxi(kk)*coordse(kk);
        dxcdxi = dxcdxi + dshdxi(kk)*coordsce(kk);
    end
    %
    J = abs(dxdxi);
    Jc = abs(dxcdxi);
    %
    dshdx = zeros(1,nNodese);
    dshdxc = zeros(1,nNodese);
    for kk = 1:nNodese
        dshdx(kk) = dshdxi(kk)/dxdxi;
        dshdxc(kk) = dshdxi(kk)/dxcdxi;
    end
    %
    % Calculate the radial coordinate
    %
    r = 0.;
    rc = 0;
    for kk = 1:nNodese
        r = r + sh(kk)*coordse(kk);
        rc = rc + sh(kk)*coordsce(kk);
    end
    %
    % Calculate the stretches
    %
    lambdar = 1.;
    for kk = 1:nNodese
        lambdar = lambdar + dshdx(kk)*ue(kk);
    end
    %
    lambdahoop = rc/r;
    detF = lambdar*lambdahoop^2;
    %
    % Compute the spatial velocity gradient
    %
    epsrdot = 0.;
    epshoopdot = 0.;
    for kk = 1:nNodese
        epsrdot = epsrdot + dshdxc(kk)*ve(kk);
        epshoopdot = epshoopdot + sh(kk)*ve(kk)/rc;
    end
    %
    % Evolution equation for the viscous hoop stretch: Explicit time integration
    if (dt==0)
        % Initial dummy step - initialize lambdav = lambdahoop and Me = 0, i.e., a fully-relaxed Maxwell element
        lambdav_new = lambdahoop;
        Me_new = 0.0;
    else
        lambdav_new = exp(dt*gdot0*((sqrt(3)*abs(Me_old)/shearmod)^(1/mp))*sign(Me_old)/(2*sqrt(3)))*lambdav_old; % Update the viscous hoop stretch (uses the exponential map)
        lambdae_bar = lambdahoop/lambdav_new/detF^(1/3); % Update the elastic distortional hoop stretch
        Me_new = 2*G1*log(lambdae_bar); % Update the hoop component of the elastic Mandel stress
    end
    %
    % Compute the Cauchy stress
    %
    %--- equilibrium branch - quadratic-I1 ---%
    %
    quad_fac = 1 + alpha*((lambdar^2 + 2*lambdahoop^2)/detF^(2/3) - 3 );
    sigmar_eq = (shearmod/detF^(5/3))*quad_fac*(lambdar^2 - (lambdar^2 + 2*lambdahoop^2)/3) + bulkmod*(detF-1);
    sigmahoop_eq = (shearmod/detF^(5/3))*quad_fac*(lambdahoop^2 - (lambdar^2 + 2*lambdahoop^2)/3) + bulkmod*(detF-1);
    %
    %--- non equilibrium branch - Maxwell element ---%
    sigmar_neq = -2*Me_new/detF;
    sigmahoop_neq = Me_new/detF;
    %
    %--- Total stress = eq stress + non_eq stress (+ bulk viscosity contribution if needed) ---%
    sigmar = sigmar_eq + sigmar_neq + bulkvisco*(epsrdot + 2*epshoopdot);
    sigmahoop = sigmahoop_eq + sigmahoop_neq + bulkvisco*(epsrdot + 2*epshoopdot);  
    %
    % Add contribution to element residual and mass matrix from
    % current integration point
    %
    for kk = 1:nNodese
        %
        %---------- Element residual vector -----------------%
        Re(kk) = Re(kk) - sigmar*w(jj)*Jc*dshdxc(kk)*(4*pi*rc^2) - 2*sigmahoop*w(jj)*Jc*sh(kk)/rc*(4*pi*rc^2);
        %
        %------------ Element mass matrix -------------------%
        for ll = 1:nNodese
            Me(kk,ll) = Me(kk,ll) + w(jj)*J*sh(kk)*sh(ll)*(4*pi*r^2);
        end
        
    end
    %
end

end

%%
% Function that the ODE Solver calls to march governing equations for the
% bubble contents forward in time
function dxdt = bubble(t,x,params,tstart,tend,Rstart,Rend)

% Extract quantities from the parameters vector
NT = params(1); % Mesh points inside the bubble
Rv = params(2); % Gas constant for vapor (J/kg-K)
Ra = params(3); % Gas constant for air (J/kg-K)
kappa = params(4); % Specific heats ratio
fom = params(5); % Mass Fourier number
chi = params(6); % Lockhart–Martinelli number
A_star = params(7); % Dimensionless A parameter
B_star = params(8); % Dimensionless B parameter (Note that A_star+B_star=1.)
Pv_star = params(9); % Dimensionless vapor saturation pressure at the far field temperature

% Extract the bubble velocity and radius at time t
ttilde = (t - tstart)/(tend - tstart); % Normalized time over the time increment (0<ttilde<1)
R = Rstart + ttilde*(Rend - Rstart); % Linear interpolation of R(t) over the time increment
U = (Rend - Rstart)/(tend-tstart); % U(t) is calculated as dR/dt

% Extract quantities from the state vector
P = x(1); % Internal bubble pressure
Theta = x(2:(NT+1)); % Variable relating to internal temp (theta)
k = x((NT+2):(2*NT+1)); % Vapor mass fraction (k)

%******************************************
% Set up grid inside the bubble
deltaY = 1/(NT-1); % Dimensionless grid spacing inside the bubble
ii = 1:1:NT;
yk = ((ii-1)*deltaY)'; % Dimensionless grid points inside the bubble
%******************************************

%******************************************
% Apply the Dirichlet BC for the vapor mass fraction at the bubble wall
k(end) = (1+(Rv/Ra)*(P/Pv_star-1))^(-1);
%******************************************

%******************************************
% Calculate mixture fields inside the bubble
T = (A_star - 1 + sqrt(1+2.*A_star.*Theta))./A_star; % Dimensionless temperature T/T_inf
K_star = A_star.*T+B_star; % Dimensionless mixture thermal conductivity field
Rmix = k.*Rv + (1-k).*Ra; % Mixture gas constant field (J/kg-K)
%******************************************

%******************************************
% Calculate spatial derivatives of the temp and vapor conc fields
DTheta = [0; % Neumann BC at origin
    (Theta(3:end)-Theta(1:end-2))/(2*deltaY); % Central difference approximation for interior points
    (3*Theta(end)-4*Theta(end-1)+Theta(end-2))/(2*deltaY)]; % Backward difference approximation at the bubble wall
DDTheta = [6*(Theta(2)-Theta(1))/deltaY^2; % Laplacian in spherical coords at the origin obtained using L'Hopital's rule
    (diff(diff(Theta)/deltaY)/deltaY + (2./yk(2:end-1)).*DTheta(2:end-1)); % Central difference approximation for Laplacian in spherical coords
    ((2*Theta(end)-5*Theta(end-1)+4*Theta(end-2)-Theta(end-3))/deltaY^2+(2/yk(end))*DTheta(end))]; % Laplacian at the bubble wall does not affect the solution
Dk = [0; % Neumann BC at origin
    (k(3:end)-k(1:end-2))/(2*deltaY); % Central difference approximation for interior points
    (3*k(end)-4*k(end-1)+k(end-2))/(2*deltaY)]; % Backward difference approximation at the bubble wall
DDk = [6*(k(2)-k(1))/deltaY^2; % Laplacian in spherical coords at the origin obtained using L'Hopital's rule
    (diff(diff(k)/deltaY)/deltaY + (2./yk(2:end-1)).*Dk(2:end-1)); % Central difference approximation for Laplacian in spherical coords
    ((2*k(end)-5*k(end-1)+4*k(end-2)-k(end-3))/deltaY^2+(2/yk(end))*Dk(end))]; % Laplacian at the bubble wall does not affect the solution
%******************************************

%******************************************
% Internal bubble pressure evolution equation
pdot = 3/R*(-kappa*P*U + (kappa-1)*chi*DTheta(end)/R ...
    + kappa*P*fom*Rv*Dk(end)/(R*Rmix(end)*(1-k(end))));
%******************************************

%******************************************
% Dimensionless mixture velocity field inside the bubble
Umix = ((kappa-1).*chi./R.*DTheta-R.*yk.*pdot./3)./(kappa.*P) + fom./R.*(Rv-Ra)./Rmix.*Dk;
%******************************************

%******************************************
% Evolution equation for the temperature (theta) of the mixture inside the bubble
Theta_prime = (pdot + (DDTheta).*chi./R.^2).*(K_star.*T./P.*(kappa-1)./kappa) ...
    - DTheta.*(Umix-yk.*U)./R ...
    + fom./(R.^2).*(Rv-Ra)./Rmix.*Dk.*DTheta;
Theta_prime(end) = 0; % Dirichlet BC at the bubble wall
%******************************************

%******************************************
% Evolution equation for the vapor concentration inside the bubble
k_prime = fom./R.^2.*(DDk + Dk.*(-((Rv - Ra)./Rmix).*Dk - DTheta./sqrt(1+2.*A_star.*Theta)./T)) ...
    - (Umix-U.*yk)./R.*Dk;
k_prime(end) = 0; % Dirichlet BC at the bubble wall
%******************************************

dxdt = [pdot; Theta_prime; k_prime];

end