%%%%%%%%%%%%%%%%%%%%%%% Prelims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%%
%%%%%%%%%%%%%%%%%%%%%%% Endogenous Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var
c		% consumption 
n		% employment 
k		% capital
y		% output
z		% TFP
phi1    % the disutility of labor

% auxiliary variables
kn		% capital stock per employee
cn		% consumption per capita
rk		% rental rate of capital
R		% gross real interest rate
w		% wage rate
i		% investment
;

%%
%%%%%%%%%%%%%%%%%%%%%%% Exogenous Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
eps	        % productivity shocks	
sig 
;  
    
%%
%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters beta, alpha, delta, sigma, phi, rho, rhop;

beta	= 0.99; 	                      % discount rate (TO DO)
alpha	= 1/3;		                      % capital share (TO DO)
delta	= 0.025;	                      % depreciation rate
sigma	= 10;		                      % inverse elasticity of inter-temporal substitution (TO DO)
phi		= 1;	 	                      % inverse Frisch elasticity of labor supply (TO DO)
rho		= 0.9;                            % persisitence of TFP shocks
rhop    = 0.5;		                 


%%
%%%%%%%%%%%%%%%%%%%%%%% Model Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

% Households
R = rk + 1 - delta;                        % Gross real rate
c^(-sigma) = beta * c(+1)^(-sigma) * R(+1);   % Euler equation (TO DO)
n^phi  = w * c^(-sigma)/phi1;                        % Labor supply (TO DO)
cn		= c / n;                       % Consumption per capita

% Firms
y		= z*(k(-1)^alpha)*(n^(1-alpha)); % Production function
w		= (1-alpha) * z *kn^alpha;       % Wage rate
rk		= alpha * kn^(alpha-1);    % Rental rate of capital (TO DO)
kn		= k(-1) / n;                   % Capital stock per employee


% Market clearing
c + k = w*n + rk*k + (1 - delta)*k(-1);           % Good market clearing (TO DO)

% Get investment from capital lom
i = k - (1-delta)*k(-1);               % Investment 

% Shocks
log(z) =  rho * log(z(-1)) + eps;                              % Technology process (TO DO)
log(phi1) = rhop*log(phi1(-1)) + sig;
end;

%%%%%%%%%%%%%%%%%%%%%%% Steady State %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;
z	= 1;
phi1 = 1;                                
R	= 1/beta;
rk	= R - (1 - delta);
kn	= (rk/alpha)^(1/(alpha-1));
w	= (1-alpha)*(kn^alpha);
cn	= kn^alpha - delta*kn;
n	= ( (cn^(-sigma))*(1-alpha)*(kn^alpha)/phi1)^( 1 / (phi+sigma) );
k	= kn*n;
c	= cn*n;
y	= (kn^alpha)*n;
i	= delta*k;
end;

resid;
steady;
check;

%%%%%%%%%%%%%%%%%%%%%%%%% Shocks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;                 
var sig = 0.01^2;    % Unit Shock to labor                        
end;

%% Obtain policy functions, IRFs and second moments
stoch_simul(irf=20,loglinear,order=1,hp_filter=1600)
y c n k w z;                         


