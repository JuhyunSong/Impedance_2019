
%% May 2019; Juhyun Song
% Copyright is granted to the publisher of J.Song, E.Khoo, and M.Z.Bazant, 2019.
% This script calculates and plots the steady state solution of the symmetric model considered in J.Song, E.Khoo, and M.Z.Bazant, 2019.
% Written in MATLAB version 2017b.
% Requires Lambert_W.m. Download at (https://www.mathworks.com/matlabcentral/fileexchange/43419-the-lambert-w-function, last accessed 05.15.2019)
% Specify input parameters in INPUT section and run the script.


clear; clc; close all;
%% INPUTS
rho = -0.01;    % dimensionless charge density of immobile charge
j = 3;          % dimensionless current (steady state)
x = linspace(0,1,100); % dimensionless position vector for calculation

%% CHECK THE VALIDITY OF THE INPUTS
% If rho >= 0, there is a limit on current.
if rho == 0
    jmax = 2;
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


%% ANALYTICAL SOLUTION

% First try the explicit form
[c,p,dcdx,dpdx,~,~] = func_steady_explicit(x,v,rho);
y = [c;p;dcdx;dpdx];

% If the explicit form does not work, fall back to the implicit form.
if sum(sum(isnan(y))) > 0 || sum(sum(isinf(y))) > 0
    fprintf('The explicit solution does not work. Fall back to the implicit solution. \n')
    option_implicit = optimset('TolX',1e-10); dx = 10^-5;
    [c,p,dcdx,dpdx,~] = func_steady_implicit(x,v,rho,dx,option_implicit);
    y = [c;p;dcdx;dpdx];
end

%% PLOT
figure(1);
title_str = {'c(x)','p(x)','dc/dx(x)','dp/dx(x)'};

for k = 1:length(j)
    for n = 1:4
    subplot(2,2,n); plot(x,y(n,:),'-'); hold on;
    ylabel(title_str{n}); xlabel('x')
    end
end


%% ATTACHED FUNCTIONS

function [vmax,jmax] = func_maxbias(rho)

vmax = fsolve(@(v)vmax_eqn(v,rho),-log(2*rho),optimset('Display','none')); % initial guess from the reservoir configuration.
a = ((1-exp(-vmax))*(1+rho)+sqrt(((1-exp(-vmax)).^2)*((1+rho)^2)-2*rho*vmax*(1-exp(-2*vmax))))...
                    /(1-exp(-2*vmax));
jmax = a*(1 - exp(-vmax)) - rho*vmax;

        function y = vmax_eqn(v,rho)
            y = 2*rho*exp(v)...
              - ((1-exp(-v))*(1+rho)+sqrt(((1-exp(-v))^2)*((1+rho).^2)-2*rho*v*(1-exp(-2*v))))...
                     /(1-exp(-2*v));
        end

end

function [v] = func_I2V(j,rho)

if j == 0
   v = 0; 
else
   v= fsolve(@(v)i2v_eqn(v,j,rho),0.1,optimset('Display','none'));
end

        function y = i2v_eqn(v,j,rho)
            a = ((1-exp(-v)).*(1+rho)+sqrt(((1-exp(-v)).^2)*((1+rho)^2)-2*rho*v.*(1-exp(-2*v))))...
                ./(1-exp(-2*v));
            y = -j + (1-exp(-v)).*a-rho*v;
        end
    
end

function [c,p,dcdx,dpdx,d2cdx2,d2pdx2] = func_steady_explicit(x,v,rho)

if v == 0
    a = 1;
else
    a = ((1-exp(-v))*(1+rho)+sqrt(((1-exp(-v)).^2)*((1+rho)^2)-2*rho*v*(1-exp(-2*v))))...
        /(1-exp(-2*v));
end

j = a*(1 - exp(-v)) - rho*v;

if rho == 0 
    c = a-j*x;
    p = log(c/a);
    dcdx = -j;
    d2cdx2 = 0;
    dpdx = -j./c;
    d2pdx2 = -(j^2)./(c.^2);
    
elseif rho < 0
    c = -rho*Lambert_W(a/(-rho)*exp((a-j*x)/(-rho)),0);
    p = log(c/a);
    dcdx = -j*c./(c-rho);
    d2cdx2 = dcdx*(rho*j)./((c-rho).^2);
    dpdx = -j./(c-rho);
    d2pdx2 = dcdx.*j./((c-rho).^2);
elseif rho > 0
    c = -rho*Lambert_W(a/(-rho)*exp((a-j*x)/(-rho)),-1);
    p = log(c/a);
    dcdx = -j*c./(c-rho);
    d2cdx2 = dcdx*(rho*j)./((c-rho).^2);
    dpdx = -j./(c-rho);
    d2pdx2 = dcdx.*j./((c-rho).^2);
end

end

function [c,p,dcdx,dpdx,d2pdx2] = func_steady_implicit(x,v,rho,dx,options)

if v == 0 
    a =1;
else
    a = ((1-exp(-v))*(1+rho)+sqrt(((1-exp(-v)).^2)*((1+rho)^2)-2*rho*v*(1-exp(-2*v))))...
        /(1-exp(-2*v));
end
j=(1-exp(-v))*a-rho*v; % current


%% Calculating solutions
M=length(x); % resolution number in space
c=zeros (1,M); % reserving a vector
p=zeros (1,M); % reserving a vector
p0=-1.01*v; % c bottom bound for the implicit solver 
p1=0; % c upper bound for the implicit solver; theoretical bound = a


if rho == 0
    for m=1:M
        c(m)=a-j*x(m);
        p(m)=log(c(m)/a);
    end
else
    for m=1:M
        p(m)=fminbnd(@(phi)potential(rho,a,j,x(m),phi),p0,p1,options);
        c(m)=a*exp(p(m));
    end
end

%% Calculating Derivatives
p_1 = zeros(1,M); p_2 = zeros(1,M); p_3 = zeros(1,M);
c_1 = zeros(1,M); c_2 = zeros(1,M); c_3 = zeros(1,M);
dpdx = zeros(1,M); dcdx = zeros(1,M); d2pdx2 = zeros(1,M);
for m = 1:M
        if x(m) <= dx        
        p_1(m) = p(m);
        p_2(m) = fminbnd(@(p)potential(rho,a,j,x(m)+dx,p),p0,p1,options);
        p_3(m) = fminbnd(@(p)potential(rho,a,j,x(m)+2*dx,p),p0,p1,options);
        c_1(m) = c(m);
        c_2(m) = a*exp(p_2(m));
        
        dpdx(m)=(p_2(m)-p_1(m))/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_2(m)-c_1(m))/dx; 
         
        elseif x(m) >= 1-dx
        p_1(m) = fminbnd(@(p)potential(rho,a,j,x(m)-2*dx,p),p0,p1,options);
        p_2(m) = fminbnd(@(p)potential(rho,a,j,x(m)-dx,p),p0,p1,options);
        p_3(m) = p(m);
        c_2(m) = a*exp(p_2(m));
        c_3(m) = c(m);  
        
        dpdx(m)=(p_3(m)-p_2(m))/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_3(m)-c_2(m))/dx;

        else
        p_1(m) = fminbnd(@(p)potential(rho,a,j,x(m)-dx,p),p0,p1,options);
        p_2(m) = p(m);
        p_3(m) = fminbnd(@(p)potential(rho,a,j,x(m)+dx,p),p0,p1,options);
        c_1(m) = a*exp(p_1(m));
        c_2(m) = c(m);
        c_3(m) = a*exp(p_3(m)); 
        
        dpdx(m)=(p_3(m)-p_1(m))/2/dx;
        d2pdx2(m)=(p_3(m)-2*p_2(m)+p_1(m))/dx^2;
        dcdx(m)=(c_3(m)-c_1(m))/2/dx;

        end
end
        function f=potential(rho_s,a,j,x,phi)
        f=(a*(1-exp(phi))-j*x+rho_s*phi)^2;
        end

end