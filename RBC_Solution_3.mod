% Solutions - Exercise 3 Part 3

%%%%%%%%%%%%%%%%%%%%%%%%% EX.3.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Endogenous Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NEW: added fi variable. 

var
c		% consumption 
n		% employment 
k		% capital
y		% output
z		% TFP
varphi      % disutility of labor

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

% Now includes shock to disutility 

varexo 
eps		% productivity shocks
st      % shock to the disutility of labor varphi
;  
    
%%
%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now includes rho_varphi


parameters beta, alpha, delta, sigma, phi, rho, rho_varphi, sigma_eta;

beta	= 0.99; 	                  % discount rate 
alpha	= 1/3;		                  % capital share 
delta	= 0.025;	                  % depreciation rate
sigma	= 1;		                  % inverse elasticity of inter-temporal substitution 
phi		= 1;	 	                  % inverse Frisch elasticity of labor supply 
rho		= 0.9;		                  % persistence of TFP shocks
rho_varphi  = 0.5;                        % persistence of disutility shock


%%
%%%%%%%%%%%%%%%%%%%%%%% Model Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Labor supply with shock + LOM for shock 


model;

% Households
R = rk + 1 - delta;                        % Gross real rate
c^(-sigma) = beta * (c(+1)^(-sigma) * R);  % Euler equation 
w * c^(-sigma) = varphi*n^(phi);               % Labor supply
cn = c / n;                                % Consumption per capita

% Firms
y		= z*(k(-1)^alpha)*(n^(1-alpha)); % Production function
w		= (1-alpha) * z *kn^alpha;       % Wage rate
rk		= alpha * z * kn^(alpha-1);      % Rental rate of capital
kn		= k(-1) / n;                     % Capital stock per employee


% Market clearing
c + k = (1-delta)*k(-1) + w*n + rk*k(-1);  % Good market clearing 

% Get investment from capital lom
i = k - (1-delta)*k(-1);               % Investment 

% Shocks
log(z) = rho * log(z(-1)) + eps;       % Technology process
log(varphi) = rho_varphi * log(varphi(-1)) + st;   % Shock process for Disutility

end;

%%%%%%%%%%%%%%%%%%%%%%% Steady State %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NEW: added varphi=1 in steady state.

steady_state_model;
z	= 1; 
varphi  = 1;                    % Disutility of labor shock in steady state                               
R	= 1/beta;
rk	= R - (1 - delta);
kn	= (rk/alpha)^(1/(alpha-1));
w	= (1-alpha)*(kn^alpha);
cn	= kn^alpha - delta*kn;
n	= ( (cn^(-sigma))*(1-alpha)*(kn^alpha) )^( 1 / (phi+sigma) );
k	= kn*n;
c	= cn*n;
y	= (kn^alpha)*n;
i	= delta*k;
end;

steady;
check;

%%%%%%%%%%%%%%%%%%%%%%%%% Shocks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Added shock to varphi. Comment out technology shock 

shocks;
%var eps  = 1^2;                   % Unit Shock to technology
var st    = 1^2;                   % Unit shock to disutility of labor
end;

%% Obtain policy functions, IRFs and second moments
stoch_simul(irf=40,loglinear,order=1,hp_filter=1600)
y c n k w z varphi;                       


%%%%%%%%%%%%%%%%%%%%%%%%% Ex.3.2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot IRFs for key variables (Shock to varphi)

figure;

% IRF of Output (y)
subplot(2, 3, 1);
plot(oo_.irfs.y_st, 'LineWidth', 2);
title('IRF of Output (y)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% IRF of Consumption (c)
subplot(2, 3, 2);
plot(oo_.irfs.c_st, 'LineWidth', 2);
title('IRF of Consumption (c)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% IRF of Labor (n)
subplot(2, 3, 3);
plot(oo_.irfs.n_st, 'LineWidth', 2);
title('IRF of Labor (n)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% IRF of Capital (k)
subplot(2, 3, 4);
plot(oo_.irfs.k_st, 'LineWidth', 2);
title('IRF of Capital (k)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% IRF of Wages (w)
subplot(2, 3, 5);
plot(oo_.irfs.w_st, 'LineWidth', 2);
title('IRF of Wages (w)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% IRF of Disutility of Labor (varphi)
subplot(2, 3, 6);
plot(oo_.irfs.varphi_st, 'LineWidth', 2);
title('IRF of Disutility of Labor (varphi)');
xlabel('Periods');
ylabel('Deviation');
grid on;

% Add a shared title
sgtitle('IRFs for 1% Positive Shock to Disutility of Labor (varphi)');