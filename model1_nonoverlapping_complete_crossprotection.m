%% Non-overlapping Commpartmental Formulation with Complete Cross-Protection
% Deterministic four-strain multilocus immune-history model

% Partial Epitope Coverage = PEC = P1 and P2, Full Epitope Coverage = FEC = P3

% Strains are defined by two loci: ax, ay, bx, by

% Immune outcomes are channeled through the following probabilities:
% P1 = locus 1 / PEC response: a or b
% P2 = locus 2 / PEC response: x or y
% P3 = locus 1 and locus 2 / FEC response: ax, ay, bx, or by

% The model explicitly tracks all possible ordered immune histories,
% giving 32 epidemiological state variables.
% Setting up ODE function
function [T,Y] = model1_nonoverlapping_complete_crossprotection(beta_ax,beta_ay,beta_bx,beta_by,sigma_ax,sigma_ay,sigma_bx,sigma_by,P1,P2,P3,mu,time)

% Pack model parameters
pars = [beta_ax,beta_ay,beta_bx,beta_by,sigma_ax,sigma_ay,sigma_bx,sigma_by,P1,P2,P3,mu];

% Initial population vectors: y(1) = susceptible hosts, y(2:5) = initial infections with ax, ay, bx, by
% remaining compartments are initially empty immune-history states.
init_pop=[0.99799, 0.000501, 0.000502, 0.000503, 0.000504, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0, ...
    0,0,0]; 

% Specify use of ODE45 numerical solver, start at zero, steps of 3/360, ends at time
[T,Y]= ode45(@equations, 0:3/360:time,init_pop, ...
    odeset('nonnegative',1:32,'AbsTol', 1e-14), pars); 

end % end main function

% Define ODE equations
function dydt = equations(~,y,pars)

% Unpack Parameters back to original names 
beta_ax   = pars(1);
beta_ay   = pars(2);
beta_bx   = pars(3);
beta_by   = pars(4);
sigma_ax  = pars(5);
sigma_ay  = pars(6);
sigma_bx  = pars(7);
sigma_by  = pars(8);
P1        = pars(9);
P2        = pars(10);
P3        = pars(11);
mu        = pars(12);

% Unpack epidemiological state variables.
% Variable names encode infection or immune history:
% Iax = currently infected with strain ax
% Ra  = recovered with immunity to epitope a
% Rx  = recovered with immunity to epitope x
% Rax = recovered with full strain-specific immunity to ax or recovered
% with partial epitope coverage against a and x independently

S       = y(1);
Iax     = y(2);
Iay     = y(3);
Ibx     = y(4);
Iby     = y(5);
Qa      = y(6);
Qb      = y(7);
Qx      = y(8);
Qy      = y(9);
Iax_y   = y(10);
Iax_b   = y(11);
Iay_x   = y(12);
Iay_b   = y(13);
Ibx_y   = y(14);
Ibx_a   = y(15);
Iby_x   = y(16);
Iby_a   = y(17);
Rax     = y(18);
Ray     = y(19);
Rbx     = y(20);
Rby     = y(21);
Rxy     = y(22);
Rab     = y(23);
Iby_ax  = y(24);
Ibx_ay  = y(25);
Iay_bx  = y(26);
Iax_by  = y(27);
Raxy    = y(28);
Rabx    = y(29);
Raby    = y(30);
Rbxy    = y(31);
Rabxy   = y(32);

% Calculate

% Total population
N = sum(y); 

% Forces of infection CHANGE BETA TO LAMBDA
lambda_ax = beta_ax * (Iax + Iax_y + Iax_b + Iax_by) / N;
lambda_ay = beta_ay * (Iay + Iay_x + Iay_b + Iay_bx) / N;
lambda_bx = beta_bx * (Ibx + Ibx_y + Ibx_a + Ibx_ay) / N;
lambda_by = beta_by * (Iby + Iby_x + Iby_a + Iby_ax) / N;

% Define ODEs
dS   = mu*N - S*(lambda_ax + lambda_ay + lambda_bx + lambda_by) - mu*S;

% primary infections
dIax = lambda_ax*S - (sigma_ax + mu)*Iax;
dIay = lambda_ay*S - (sigma_ay + mu)*Iay;
dIbx = lambda_bx*S - (sigma_bx + mu)*Ibx;
dIby = lambda_by*S - (sigma_by + mu)*Iby;


% quasi-immune single-epitope states
dQa  = sigma_ax*Iax*P1 + sigma_ay*Iay*P1 - (lambda_bx*Qa + lambda_by*Qa + mu*Qa);
dQb  = sigma_bx*Ibx*P1 + sigma_by*Iby*P1 - (lambda_ax*Qb + lambda_ay*Qb + mu*Qb);
dQx  = sigma_ax*Iax*P2 + sigma_bx*Ibx*P2 - (lambda_ay*Qx + lambda_by*Qx + mu*Qx); 
dQy  = sigma_ay*Iay*P2 + sigma_by*Iby*P2 - (lambda_ax*Qy + lambda_bx*Qy + mu*Qy); 

% secondary infections from non-blocking Q
dIax_y = lambda_ax*Qy - (sigma_ax + mu)*Iax_y;
dIax_b = lambda_ax*Qb - (sigma_ax + mu)*Iax_b;
dIay_x = lambda_ay*Qx - (sigma_ay + mu)*Iay_x;
dIay_b = lambda_ay*Qb - (sigma_ay + mu)*Iay_b;

dIbx_y = lambda_bx*Qy - (sigma_bx + mu)*Ibx_y;
dIbx_a = lambda_bx*Qa - (sigma_bx + mu)*Ibx_a;
dIby_x = lambda_by*Qx - (sigma_by + mu)*Iby_x;
dIby_a = lambda_by*Qa - (sigma_by + mu)*Iby_a;

% dual immunity 
dRax = sigma_ax*Iax*P3 + sigma_ay*Iay_x*P1 + sigma_bx*Ibx_a*P2 - (lambda_by*Rax + mu*Rax);
dRay = sigma_ay*Iay*P3 + sigma_ax*Iax_y*P1 + sigma_by*Iby_a*P2 - (lambda_bx*Ray + mu*Ray);
dRbx = sigma_bx*Ibx*P3 + sigma_by*Iby_x*P1 + sigma_ax*Iax_b*P2 - (lambda_ay*Rbx + mu*Rbx);
dRby = sigma_by*Iby*P3 + sigma_bx*Ibx_y*P1 + sigma_ay*Iay_b*P2 - (lambda_ax*Rby + mu*Rby);
dRxy = sigma_ax*Iax_y*P2 + sigma_ay*Iay_x*P2 + sigma_bx*Ibx_y*P2 + sigma_by*Iby_x*P2 - mu*Rxy;
dRab = sigma_ax*Iax_b*P1 + sigma_ay*Iay_b*P1 + sigma_bx*Ibx_a*P1 + sigma_by*Iby_a*P1 - mu*Rab;

dIby_ax = lambda_by*Rax - (sigma_by + mu)*Iby_ax;
dIbx_ay = lambda_bx*Ray - (sigma_bx + mu)*Ibx_ay;
dIay_bx = lambda_ay*Rbx - (sigma_ay + mu)*Iay_bx;
dIax_by = lambda_ax*Rby - (sigma_ax + mu)*Iax_by;

% triple immunity
dRaxy = sigma_ay*Iay_x*P3 + sigma_ax*Iax_y*P3 + sigma_by*Iby_ax*P2 + sigma_bx*Ibx_ay*P2 - mu*Raxy;
dRabx = sigma_bx*Ibx_a*P3 + sigma_by*Iby_ax*P1 + sigma_ax*Iax_b*P3 + sigma_ay*Iay_bx*P1 - mu*Rabx;
dRaby = sigma_ay*Iay_b*P3 + sigma_bx*Ibx_ay*P1 + sigma_by*Iby_a*P3 + sigma_ax*Iax_by*P1 - mu*Raby;
dRbxy = sigma_by*Iby_x*P3 + sigma_bx*Ibx_y*P3 + sigma_ay*Iay_bx*P2 + sigma_ax*Iax_by*P2 - mu*Rbxy;

dRabxy = sigma_by*Iby_ax*P3 + sigma_bx*Ibx_ay*P3 + sigma_ay*Iay_bx*P3 + sigma_ax*Iax_by*P3 - mu*Rabxy;


% pack derivatives
dydt = [dS; dIax; dIay; dIbx; dIby; ...
        dQa; dQb; dQx; dQy; ...
        dIax_y; dIax_b; dIay_x; dIay_b; ...
        dIbx_y; dIbx_a; dIby_x; dIby_a; ...
        dRax; dRay; dRbx; dRby; dRxy; dRab; ...
        dIby_ax; dIbx_ay; dIay_bx; dIax_by; ...
        dRaxy; dRabx; dRaby; dRbxy; dRabxy];
end