% ResistorBallistic.m 
% Calculation of ballistic electron transport in n++-n+-n++ GaAs 
% resistor using NEGF formalism (includes transverse momentum)
% Description
% 1D GaAs 142nm long resistor structure n++ - n+ - n++ is considered. 
% n++ side sections modeled with 15 grid points, the central part consists
% of 90 grid points. Bias value is Vp = 0.2V. The code solves NEGF
% equations and Shrodinger equation self-consistently. Calculates
% potential, LDOS, current, electron density and quasi-Fermi levels.
%
%
% Syntax:  ResistorBallistic
%
% Outputs:
%    output1 - 1D plot of the potential U(x)
%    output2 - 2D plot of the LDOS(x,E)
%    output3 - 1D plot of the total electron density n(x) and doping
%              profile Nd(x)
%    output4 - plot of the spectral electron density n(x,E)
%    output5 - plot of the spectral current I(x,E)
%    output6 - plot of the quasi-Fermi levels vs x
% 
%
% Other m-files required: FermiDiracInt.m 
%                         calculates Fermi-Dirac distribution
% Subfunctions: none
% MAT-files required: 
% Author: Vladimir Gerasik,
% Wilfrid Laurier University, Dept. of Physics & CS,
% December 2014; Last revision: 29-April-2017
clear all


% Fundamental Physical Constants
h    =   6.582119514e-16;    % Plank's constant [eV*s]
hh   =   1.0545718e-34;      % Plank's constant [J*s]
m0   =   9.11e-31;           % electron mass
q    =   1.6e-19;            % elementary charge [C]
kT   =   8.6173e-5*300;      % Boltzmann constant x temperature T=300K [eV]

% Parameters GaAs
m_e   = 0.067*m0;            % effective mass [kg]
eps   = 10.5*8.85e-12;       % permittivity [F/m]

% Structure (resistor) 
Ns     =  15;               % side section 
Nc     =  90;               % central part 
Nx     =  Ns + Nc + Ns;     % total number of discrete points  
dx     =  1.2*1e-9;         % delta x [m] (device lenght is dx*Nx)

% Energy discretization
E_C = 1.420;                % band gap energy GaAs
Emin=E_C-0.22;              % energy interval 
Emax=E_C+0.58;

NE = 300;                   % number of energy grid points
E = linspace(Emin, Emax, NE); 
DE=(E(2)-E(1));             % energy spacing  

% Fermi level (0.4eV above the bandedge)
mu = E_C + 0.4; 

% Doping n++ n+ n++
NC     =  2*((m_e*kT*q/(2*pi*(hh^2)))^1.5);            % effective density of states in the conduction band
nd     =  NC*FermiDiracInt(1/2,(mu-E_C)/kT);           % reference doping level (outer regions)
Nd     =  nd*[ones(Ns,1);.1*ones(Nc,1);ones(Ns,1)]';   % doping profile of the structure (discretized)
    

% Tridiagonal matrix for Poisson equation
   D2 = -2*diag(ones(1,Nx))+...
                        diag(ones(1,Nx-1),1)+...
                               diag(ones(1,Nx-1),-1);
% Boundary conditions (Neumann homogeneous)
   D2(1,2)=2; D2(Nx,Nx-1)=2;   
   %D2(1,1)=-1; D2(Nx,Nx)=-1;  %zero field condition (Datta)

% k-points discretization
k_max = (1.5/1e-9);       % maximum k-value
NK = 80;                  % number of k-points
k_point = linspace(0,k_max,NK);
dk = k_point(2) - k_point(1);

% Transverse energies E_k = h^2k^2/2m0
E_k = h*hh*k_point.^2./(2*m_e);

% H and W matrices (bulk)
t = h*hh/(2*m_e*dx^2);   % t parameter

H0 = zeros(1,NK);

H0(1:NK) =  E_C + E_k(1:NK) + 2*t;  
W        =  -t;           
    
% Construct the main Hamiltonian matrix T
T = zeros(Nx,Nx,NK);

for k=1:NK
T(:,:,k) = H0(k)*diag(ones(1,Nx))+...
                 W*diag(ones(1,Nx-1),1)+...
                       W*diag(ones(1,Nx-1),-1);
end
%Gamma for calculations
    gamm = 1i*1e-15;    % complex parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bias
Vp = 0.20;      % bias [V]


% Fermi distrubution 
     
    f0L(1:NE) = 1./(1+exp((E(1:NE) - mu  )./kT)); 
    f0R(1:NE) = 1./(1+exp((E(1:NE) - mu + Vp)./kT));

 
% Initial potential profile (zero)
   U  =  (0+0*linspace(0,Nx,Nx)/Nx)';
   load 'Upot.mat' U

  
% Initialize Green's functions
   
  G_R         = zeros(Nx,Nx,NK);    % retarded GF
  G_L         = zeros(Nx,Nx,NK);    % lesser (electron) GF
  BSE_L       = zeros(Nx,Nx,NK);    % boundary self-energy left
  BSE_R       = zeros(Nx,Nx,NK);    % boundary self-energy right
  BSE_lesser  = zeros(Nx,Nx,NK);    % boundary self-energy lesser
  
  
  GamL = zeros(Nx,Nx,NK);           % level broadening (left)
  GamR = zeros(Nx,Nx,NK);           % level broadening (right)
  
  
  
  LDOS_en =  zeros(Nx,NE);           % LDOS
  ElDen   =  zeros(Nx,NE);           % spectral electron density
  
% Initial tolerance and outer loop step count
  TOL  = 100;
  Step = 1;
 
% MAIN CYCLE (outer cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic  
 

  while TOL>0.01
 
      
  J_C     =  zeros(Nx-1,NE);
  n       =  zeros(1,Nx);
      
    
% ENERGY CYCLE (inner cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for EnergyNode=1:NE
    
  fprintf('Main Cycle %d Energy cycle %d and convergence %d  \n',Step, EnergyNode, TOL);
   
 % Analytic calculation of k_L, k_R
 k_L(1:NK) = (1/(dx)).*acos(1 - m_e*dx^2*(E(EnergyNode) - E_C - E_k(1:NK) - U(1))/(h*hh)); 
 k_R(1:NK) = (1/(dx)).*acos(1 - m_e*dx^2*(E(EnergyNode) - E_C - E_k(1:NK) - U(Nx))/(h*hh));
 
 sig_L(1:NK) = - t*exp(1i*k_L(1:NK)*dx);
 sig_R(1:NK) = - t*exp(1i*k_R(1:NK)*dx);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Boundary self-energies 
 BSE_L(1,1,1:NK)    = sig_L(1:NK);
 BSE_R(Nx,Nx,1:NK)  = sig_R(1:NK);
 
 % Level broadening 
 for k=1:NK
    GamL(:,:,k) = 1i*(BSE_L(:,:,k) - BSE_L(:,:,k)');
    GamR(:,:,k) = 1i*(BSE_R(:,:,k) - BSE_R(:,:,k)');
 end
 
    
 % Sigma-lesser Green's function
  for k=1:NK
    BSE_lesser(:,:,k)    =   1i * (  GamL(:,:,k)*f0L(EnergyNode)  +  GamR(:,:,k)*f0R(EnergyNode) );            
  end
 
 
 % Ballistic Green's functions integrated wrt k
 G_L_T       =  zeros(Nx,Nx);  
 G_R_T       =  zeros(Nx,Nx); 
 
 % Green's function calculation for all k-points
 for k=1:NK

  G_R(:,:,k) = inv( (E(EnergyNode)+gamm)*eye(Nx) - T(:,:,k) - diag(U) - BSE_L(:,:,k) - BSE_R(:,:,k) );
  G_L(:,:,k) = G_R(:,:,k)*BSE_lesser(:,:,k)*G_R(:,:,k)';
  G_L_T = G_L_T + 1/pi*dk*k_point(k)*G_L(:,:,k);
  G_R_T = G_R_T + 1/pi*dk*k_point(k)*G_R(:,:,k);
  
 end
 
 % Electron density (total)
 
  n  = n - real(1i*DE/(2*pi*dx)*diag(G_L_T))'; 
 

 % Spectral electron density 
 
  ElDen(:,EnergyNode) = - real(1i/dx*diag(G_L_T))';
  
             
 % Local density of states 
 
         for SpaceGrid=1:Nx
                 LDOS_en(SpaceGrid,EnergyNode) = - 1/(pi*dx) * imag(G_R_T(SpaceGrid,SpaceGrid));
         end    
         
 % Spectral current
         for SpaceGrid=1:Nx-1
                 J_C(SpaceGrid, EnergyNode) = J_C(SpaceGrid,EnergyNode) + (q/h)*(W*G_L_T(SpaceGrid+1,SpaceGrid) - ... 
                           W*G_L_T(SpaceGrid,SpaceGrid+1));                     
         end
         
  
 
end

% END ENERGY CYCLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quazi-Fermi level Eq.(37)
    LHS = (n/NC)';
    QFL = kT*( log(LHS) ./ (1 - LHS.^2) + (3*sqrt(pi).*LHS./4).^(2/3) ./...
        (1 + (0.24 + 1.08*(3*sqrt(pi).*LHS/4 ).^(2/3)).^(-2) ))' + U';
    
% Averaging boundary points
    QFL(Ns)=(QFL(Ns-1)+QFL(Ns+1))/2;
    QFL(Ns+Nc+1)=(QFL(Ns+Nc)+QFL(Ns+Nc+2))/2;
            
            Step = Step+1;
   
 % Newton-Raphson iteration (correction dU from Poisson equation) 
             
           
            D=zeros(Nx,1); 
            for k=1:Nx 
                        z=(QFL(k) - U(k))/kT; 
                        D(k,1)=q*dx^2/(eps*kT)*NC*FermiDiracInt(-1/2,z); 
            end 
                        dN=-(q*dx^2/eps)*(n-Nd)-(D2*U)'; 
                        dU=(D2-diag(D))\dN';
                        
                         U=U+dU; 
                              
                        TOL = norm(dU)/(norm(U));
                        
             
end
% END OUTER CYCLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\nFinal convergence %0.4f  \n', TOL);



toc


% Potential plot
 figure('Name','Electrostatic potential','NumberTitle','off'); 
                 plot(0.5:1:Nx-0.5,U+E_C)


% Electron Density plot
  figure('Name','Electron Density','NumberTitle','off');          
          plot(n)
            hold on
               plot(Nd)   
                  hold off
            
            
            
            
ntick=8;
% Electron Spectral Density plot           
figure('Name','Spectral Electron Density','NumberTitle','off');
colormap(hot);
clims = [0 4e+26];
imagesc(ElDen(:,:)')
set(gca,'YDir','normal');
set(gca,'YTick',1:(NE-1)/ntick:NE);
set(gca,'YTickLabel',E(1):(E(NE)-E(1))/ntick:E(NE) ); 
ylabel('E[eV]');          
      


% LDOS plot
figure('Name','Local Density of States','NumberTitle','off'); 
colormap(hot);
imagesc(LDOS_en(:,:)')
set(gca,'YDir','normal');
set(gca,'YTick',1:(NE-1)/ntick:NE);
set(gca,'YTickLabel',E(1):(E(NE)-E(1))/ntick:E(NE)); 
ylabel('E[eV]'); 

% Current plot
figure('Name','Spectral Current','NumberTitle','off'); 
colormap(hot);
imagesc(real(J_C(:,:)'))
numpxl=size(J_C);
set(gca,'YDir','normal');
set(gca,'YTick',1:(NE-1)/ntick:NE);
set(gca,'YTickLabel',E(1):(E(NE)-E(1))/ntick:E(NE)); 
set(gca,'XTick',1:(numpxl(1)-1)/(ntick-2):numpxl(1));
set(gca,'XTickLabel',0:(Nx)/(ntick-2):Nx);
ylabel('E[eV]'); 

J_C_TOTAL=0;
for e=1:NE
    J_C_TOTAL = J_C_TOTAL +  DE/(2*pi)*real(J_C(20,e));
end


% 3D-DOS at point x plot vs analytic 
x=20;
figure('Name','LDOS calculated vs analytic','NumberTitle','off');
plot(E(1:NE),LDOS_en(x,1:NE)');
hold on
plot(E(1:NE), (2*m_e/(hh*h))^(3/2)/(2*pi^2)*real(sqrt(E(1:NE) - E_C - U(x))));

% Plot of quasi-Fermi level
figure('Name','Quasi Fermi level','NumberTitle','off'); 
plot(0.5:1:Nx-0.5,QFL)


