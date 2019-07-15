
%% May 2019; Juhyun Song
% Copyright is granted to the publisher of J.Song, E.Khoo, and M.Z.Bazant, 2019.
% This script calculates and plots the impedanec solution of the reservoir model considered in J.Song, E.Khoo, and M.Z.Bazant, 2019.
% Written in MATLAB version 2017b.
% Requires Lambert_W.m. Download at (https://www.mathworks.com/matlabcentral/fileexchange/43419-the-lambert-w-function, last accessed 05.15.2019)
% Specify input parameters in INPUT section and run the script.

clear; clc; close all;
%% INPUTS
rho = -0.01;        % dimensionless charge density of immobile charge
j = 1.5;            % dimensionless current (steady state)
% frequency range
log_w_max = 7;      % log10 of the maximum frequency (Hz)
log_w_min = -2;     % log10 of the maximum frequency (Hz)
ppd = 5;            % number of points per decade 


%% CHECK THE VALIDITY OF THE INPUTS
% If rho >= 0, there is a limit on current.
if rho == 0
    jmax = 1;
elseif rho > 0
    [~,jmax] = func_maxbias(rho);
elseif rho < 0
    jmax = inf;
end
% Check if the given curent is out of the range.
if j >= jmax
    fprintf('The current is out of the valid range, 0 <= j < %2.1f .\n',jmax)
    return 
end

%% CALCULATE THE VOLTAGE
v = func_I2V(j,rho);
fprintf('For the given rho (%2.3f) and current (%2.2f), the corresponding voltage is %2.2f . \n',rho,j,v)


%% STEADY STATE SOLUTION 
% Use a numerical method to obtain the steady state solution in a specific format.
x0 = linspace(0,1,1001);
options_ss = bvpset('Vectorized','on','RelTol',1e-6,'AbsTol',1e-8,'NMax',5000);
[sol0,~] = func_steady_numerical(x0,v,rho,options_ss);

%% IMPEDANCE SOLUTION
w_exp = linspace(log_w_min,log_w_max,(log_w_max-log_w_min)*ppd+1); 
w = 10.^w_exp;
v1 = 1; % perturbation amplitude, the impedance solution should be independet of this.
options_imp = bvpset('Vectorized','on','RelTol',1e-6,'AbsTol',1e-8,'NMax',5000);
[z] = func_impedance(w,v1,sol0,rho,options_imp);  % impedance for entire w range.
[z1] = func_impedance(1,v1,sol0,rho,options_imp); % impedance at w = 1


%% PLOT
figure (1);
plot(real(z),-imag(z),'-'); hold on;
plot(real(z1),-imag(z1),'o')
% axis setting
axiscut = max(max(real(z)),max(-imag(z)));
axis([0 axiscut 0 axiscut])
% make it squared
daspect([1 1 1])



%% ATTACHED FUNCTIONS

function [vmax,jmax] = func_maxbias(rho)

vmax = -log(2*rho);
jmax = 1 - exp(-vmax) - rho*vmax;

end


function [v] = func_I2V(j,rho)

if j == 0
   v = 0; % if current is zero, voltage is zero anyhow 
else
    if rho == 0
        c1 = 1-j;
        p1 = log(c1);
        v = -p1;
    elseif rho < 0
        c1 = -rho*Lambert_W(1/(-rho).*exp((1-j)/(-rho)),0);
        p1 = log(c1);
        v = -p1;

    elseif rho > 0
        c1 = -rho*Lambert_W(1/(-rho).*exp((1-j)/(-rho)),-1);
        p1 = log(c1);
        v = -p1;
    end
end

if isinf(v) || isnan(v)
    disp('!!! alternative approach in converting j to v. \n')
    v =fsolve(@(s)func_find_v(s,j,rho),(-rho)^-1*(j-1));
end

        function [f_find_V] = func_find_v(v,j,rho)
                f_find_V = j - (1 - exp(-v) - rho*v);
        end

end


function [sol,sol_a] = func_steady_numerical(x0,v,rho,option)

sol_a = bvpinit(x0,@(x)init_guess1(x,v,rho)); % initial guess by analytical solution
    if sum(sum(isnan(sol_a.y))) > 0 || sum(sum(isinf(sol_a.y))) > 0
        disp('!!! alternative approach in initial guess of steady state numerical solution. \n')
        x0 = linspace(0,1,1001); dx = 1e-5; option_implicit = optimset('TolX',1e-10); 
        sol_a = bvpinit(x0,@(x)init_guess2(x,v,rho,dx,option_implicit));
    end
sol = bvp4c(@(x,y)gov_eqn(x,y,rho),@(y0,y1)bcd_eqn(y0,y1,v),sol_a,option);


        function y = init_guess1(x,v,rho)
        [y(1),y(2),y(3),y(4),~,~] = func_steady_explicit(x,v,rho);
        end


        function y = init_guess2(x,v,rho,dx,option_implicit)
        [y(1),y(2),y(3),y(4),~] = func_steady_implicit(x,v,rho,dx,option_implicit);
        end

        function dydx = gov_eqn(x,y,rho)
        dydx = [y(3,:)                              % c
                y(4,:)                              % p
                -rho./(y(1,:)-rho).*y(4,:).*y(3,:); % dcdx
                -1./(y(1,:)-rho).*y(3,:).*y(4,:)];  % dpdx
        end

        function y_bc = bcd_eqn(y0,y1,v)
        y_bc = [y0(1,:) - 1                         % c0 = 1
                y0(2,:)                             % p0 = 0
                y1(2,:) + v                         % p1 = -v
                y1(3,:)-y1(1,:).*y1(4,:)];          % 0 = dcdx - c*dpdx at x=1
        end

end


function [c,p,dcdx,dpdx,d2cdx2,d2pdx2] = func_steady_explicit(x,v,rho)

j = 1 - exp(-v) - rho*v;

if rho == 0 
    c = 1-j*x;
    p = log(c);
    dcdx = -j;
    d2cdx2 = 0;
    dpdx = -j./c;
    d2pdx2 = -(j^2)./(c.^2);
    
elseif rho < 0
    c = -rho*Lambert_W(1/(-rho)*exp((1-j*x)/(-rho)),0);
    p = log(c);
    dcdx = -j*c./(c-rho);
    d2cdx2 = dcdx*(rho*j)./((c-rho).^2);
    dpdx = -j./(c-rho);
    d2pdx2 = dcdx.*j./((c-rho).^2);
elseif rho > 0
    c = -rho*Lambert_W(1/(-rho)*exp((1-j*x)/(-rho)),-1);
    p = log(c);
    dcdx = -j*c./(c-rho);
    d2cdx2 = dcdx*(rho*j)./((c-rho).^2);
    dpdx = -j./(c-rho);
    d2pdx2 = dcdx.*j./((c-rho).^2);
end

end

function [c,p,dcdx,dpdx,d2pdx2] = func_steady_implicit(x,v,rho,dx,options)

j=1-exp(-v)-rho*v; % current
M = length(x); % resolution number in space


%% Calculating solutions
c=zeros (1,M); % reserving a vector
p=zeros (1,M); % reserving a vector
p0=-v; % c bottom bound for the implicit solver 
p1=0; % c upper bound for the implicit solver


if rho == 0
    for m=1:M
        c(m)=1-j*x(m);
        p(m)=log(c(m));
    end
else
    for m=1:M
        p(m)=fminbnd(@(p)potential(rho,j,x(m),p),p0,p1,options);
        c(m)=exp(p(m));        
    end
end

%% Calculating Derivatives
p_1 = zeros(1,M); p_2 = zeros(1,M); p_3 = zeros(1,M);
c_1 = zeros(1,M); c_2 = zeros(1,M); c_3 = zeros(1,M);
dpdx = zeros(1,M); dcdx = zeros(1,M); d2pdx2 = zeros(1,M);
for m = 1:M
        if x(m) <= dx        
        p_1(m) = p(m);
        p_2(m) = fminbnd(@(p)potential(rho,j,x(m)+dx,p),p0,p1,options);
        p_3(m) = fminbnd(@(p)potential(rho,j,x(m)+2*dx,p),p0,p1,options);
        c_1(m) = c(m);
        c_2(m) = exp(p_2(m));
        
        dpdx(m)=(p_2(m)-p_1(m))/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_2(m)-c_1(m))/dx;    
         
        elseif x(m) >= 1-dx
        p_1(m) = fminbnd(@(p)potential(rho,j,x(m)-2*dx,p),p0,p1,options);
        p_2(m) = fminbnd(@(p)potential(rho,j,x(m)-dx,p),p0,p1,options);
        p_3(m) = p(m);
        c_2(m) = exp(p_2(m));
        c_3(m) = c(m);  
        
        dpdx(m)=(p_3(m)-p_2(m))/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_3(m)-c_2(m))/dx;

        else
        p_1(m) = fminbnd(@(p)potential(rho,j,x(m)-dx,p),p0,p1,options);
        p_2(m) = p(m);
        p_3(m) = fminbnd(@(p)potential(rho,j,x(m)+dx,p),p0,p1,options);
        c_1(m) = exp(p_1(m));
        c_2(m) = c(m);
        c_3(m) = exp(p_3(m)); 
        
        dpdx(m)=(p_3(m)-p_1(m))/2/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_3(m)-c_1(m))/2/dx;

        end
end
        function f=potential(rho_s,j,x,p)
        f=(1-exp(p)-j*x+rho_s*p)^2;
        end
end

function [z,sol1] = func_impedance(w,v1,sol0,rho,options_imp)

%% Solving impedance for each freqeuncy
K = length(w); % number of frequency points
z = NaN(1,K); 

for k = 1:K
    % initial guess
   if k == 1 % for the first entry, use a uniform zero guess
      x_init = linspace(0,1,101);
      y_init(1,:) = zeros(size(x_init));    % c
      y_init(2,:) = v1*x_init;              % p
      y_init(3,:) = zeros(size(x_init));    % dcdx
      y_init(4,:) = v1*ones(size(x_init));  % dpdx
      init1.x = x_init;
      init1.y = y_init;
   else % then use the previous solution to be the initial guess
      init1.x = sol1.x;
      init1.y = sol1.y;      
   end
    
    % solve
      sol1 = bvp4c(@(x,y)gov_eqn1(x,y,w(k),rho,sol0),...
                   @(y0,y1)bcd_eqn1(y0,y1,v1,rho,sol0),...
                   init1,options_imp);
               
    % calculate impedance  
      y0 = deval(sol0,1);
      y1 = deval(sol1,1);
      j = (y0(1) - rho).*y1(4) + y1(1).*y0(4);
      z(k) = v1./j;
    
end
        function dydx = gov_eqn1(x,y,w,rho,sol0)
            [y_s,dy_s] = deval(sol0,x);
            dydx = [y(3,:)                              % c
                    y(4,:)                              % p
                    1i*w*y(1,:)-( y(4,:).*y_s(3,:) + y(3,:).*y_s(4,:) + y(1,:).*dy_s(4,:))*rho./(y_s(1,:)-rho) % dcdx
                    -( y(4,:).*y_s(3,:) + y(3,:).*y_s(4,:)+y(1,:).*dy_s(4,:) )./(y_s(1,:)-rho)];  % dpdx        
        end


        function y_bc = bcd_eqn1(y0,y1,v1,rho,sol0)
            y_s1 = deval(sol0,1);    
            y_bc = [y0(1,:)                             % c0 = 0
                    y0(2,:)                             % p0 = 0
                    y1(2,:) - v1                        % p1 = v1
                    y1(1,:)*y_s1(4,:) + y1(4,:)*y_s1(1,:) - y1(3,:)];          % 0 = dcdx - c*dpdx at x=1
        end
end