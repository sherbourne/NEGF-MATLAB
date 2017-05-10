% LDOS_NEGF.m - Calculation of Local Density of States (LDOS) using NEGF 
% formalism.
%
% The density of states (DOS) of a system describes the number of
% states per interval of energy at each energy level that are available 
% to be occupied. The code uses GaAs 1D nanowire parameters as an example 
% and approximates well-known analytical values of 1D LDOS and 3D LDOS 
% using NEGF formalism via retaded Green's functions.
%
% Syntax:  LDOS_NEGF
%
% Outputs:
%    output1 - 1D LDOS plot analytic and that approximated via retaded 
% Green's function
%    output2 - 3D LDOS plot analytic and that approximated via retaded 
% Green's function
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Vladimir Gerasik,
% Wilfrid Laurier University, Dept. of Physics & CS,
% December 2014; Last revision: 29-April-2017


%------------- BEGIN CODE --------------

% Physical Constants
h    =   6.582119514e-16;    % Plank's constant [eV*s]
hh   =   1.0545718e-34;      % Plank's constant [J*s]
m0   =   9.11e-31;           % Electron Mass [kg]
q    =   1.6e-19;            % Elementary Charge [C]
kT   =   8.6173e-5*300;      % Boltzmann constant kB multiplied by 
                             %           temperature T=300K

% Parameters GaAs
Eg    = 1.420;               % conduction band gap [eV]
m_eff = 0.067*m0;            % effective mass [kg]
eps   = 13.1*8.85e-12;       % permittivity [F/m]

% Structure 
Nx       =  100;             % number of spatial grid points  
dx       =  1.2*1e-9;        % delta x in nm (total device lenght is dx*Nx)


% Energy discretized

Emin = Eg-0.22;               % energy interval [eV]
Emax = Eg+0.58;

NE = 300;                      % number of energy grid points
E  = linspace(Emin, Emax, NE); % energy grid
 
gamma = 1i*1e-15;              % small complex parameter gamma Eq.(43)

% k-points (transverse momentum)
k_max = (1.5/1e-9);
NK = 100;                      % number of k grid points
k_point = linspace(0,k_max,NK);
dk = k_point(2) - k_point(1);


% transverse energies E_k = h^2k^2/2m0 Eq.(7)
E_k = h*hh*k_point.^2./(2*m_eff);

% H and W matrices (bulk Hamiltonian)
H0 = zeros(1,NK);

H0(1:NK) =  Eg + 0.5*(h*hh/m_eff)*( 2/dx^2 + k_point(1:NK).^2 );  
W        =  -0.5*(h*hh/m_eff)/dx^2;           
                       

% Construct the main Hamiltonian matrix T
T = zeros(Nx,Nx,NK);
      

for k=1:NK
      T(:,:,k) = H0(k)*diag(ones(1,Nx))+...
                        W*diag(ones(1,Nx-1),1)+...
                            W*diag(ones(1,Nx-1),-1);
end

  

  % Green's functions
   
  G_R    = zeros(Nx,Nx,NK);   % retarded Green's function
  BSE_L  = zeros(Nx,Nx,NK);   % boundary self-energy lesser (left)
  BSE_R  = zeros(Nx,Nx,NK);   % boundary self-energy lesser (right)
  
  % Level broadening functions Eq.(47)
  GamL = zeros(Nx,Nx,NK);
  GamR = zeros(Nx,Nx,NK);
  
  LDOS_3D =  zeros(Nx,NE);
  LDOS_1D =  zeros(Nx,NE);
 
tic  
% Energy cycle
fprintf('\n Energy cycle %d out of %d \n', 0, NE); 
for EnergyIter=1:NE
    
    %fprintf('Energy cycle %d out of %d \n', EnergyIter, NE); 
    
   progress = sprintf('Energy cycle %d out of %d \n', EnergyIter, NE); 
   reverseStr = repmat(sprintf('\b'), 1, length(progress));
   fprintf([reverseStr, progress]); 
    
 % Wavenumbers from dispersion relation according to Eq.(49), Eq.(48)
 kLR(1:NK) = (1/(dx)).*acos(1 - m_eff*dx^2*(E(EnergyIter) - Eg - E_k(1:NK))/(h*hh)); 
  
 
 % Boundary conditions (corner elements)
 sig_LR(1:NK) = - h*hh/(2*m_eff*dx^2)*exp(1i*kLR(1:NK)*dx);
 
 
 
 
 % Boundary self-energies lesser Eq.(48)   
 BSE_L(1,1,1:NK)                  = sig_LR(1:NK);
 BSE_R(Nx,Nx,1:NK)                = sig_LR(1:NK);
 
 
 
 % Ballistic Green's function integrated wrt k
  G_R_T       =  zeros(Nx,Nx); 
 
 for k=1:NK
  G_R(:,:,k) = inv( (E(EnergyIter)+gamma)*eye(Nx) - T(:,:,k)  - ...
                        BSE_L(:,:,k) - BSE_R(:,:,k));  % Dyson equation
  G_R_T = G_R_T + (dk/pi)*k_point(k)*G_R(:,:,k);       % G_R integrated wrt k Eq.(27)                                              
 end
 
 % LDOS approximated by Green's functions Eq.(23)
  
  for SpaceGrid=1:Nx
     LDOS_3D(SpaceGrid,EnergyIter) = - 1/(pi*dx) * imag(G_R_T(SpaceGrid,SpaceGrid));
     LDOS_1D(SpaceGrid,EnergyIter) = - 1/(pi*dx) * imag(G_R(SpaceGrid,SpaceGrid,1));
  end    
   
 
end

% END ENERGY CYCLE

toc

% 1D-DOS plot
figure('Name','1D DOS','NumberTitle','off'); 
  plot(E,LDOS_1D(2,1:NE)','DisplayName','NEGF approximation');
    hold on;
  plot(E, 1/(pi*h*sqrt(q))*real(sqrt(m_eff./(2*(E-Eg)))),...
                                                'DisplayName','Analytic expression')
       xlabel('E[eV]'); 
       ylabel('1D LDOS');
   hold off; 
       legend('show')

% 3D-DOS plot
figure('Name','3D DOS','NumberTitle','off');
  plot(E,LDOS_3D(2,1:NE)','DisplayName','NEGF approximation');
      hold on
  plot(E(1:NE), 1/(2*pi^2*(q)^(3/2))*(2*m_eff/h^2)^(3/2).*real(sqrt(E - Eg)),...
                                                'DisplayName','Analytic expression');
        xlabel('E[eV]'); 
        ylabel('3D LDOS');
      hold off;
          legend('show')
      
 %------------- END OF CODE --------------
