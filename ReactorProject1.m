clc;
close all;
clear;

% Given concentration
Ca_im = [2 0.78 0.45 0.34 0.18 0.14 0.12 0.09 0.08 0.032 0.024];
Cb_im = [1 0.78 0.64 0.50 0.36 0.28 0.22 0.17 0.14 0.09  0.048];

t = [0 1 2 3 5 7 10 15 20 40 60];   % given time

% model to calculate concentration at any given time
Ca_ic = @(alpha1) ...
        (((Ca_im(1)).^(1-alpha1(1))...
      - (1 - alpha1(1)).*alpha1(2).*t).^(1./(1-alpha1(1))));

% objective function to minimize the difference between model 
% and experimental data  
obj1 = @(alpha1) (sum((Ca_im - Ca_ic(alpha1)).^2));

alpha0 = [2 2]; % initial guess

% minimize using fmincon
opts = optimoptions(@fmincon, 'Display', 'off');
alpha1 = fmincon(obj1, alpha0, [], [], [], [], 0, [], '', opts);

% model to calculate concentration at any given time
Cb_ic = @(alpha2) ...
        (((Cb_im(1)).^(1-alpha2(1))...
      - (1 - alpha2(1)).*alpha2(2).*t).^(1./(1-alpha2(1))));

% objective function to minimize the difference between model 
% and experimental data  
obj2 = @(alpha2) (sum((Cb_im - Cb_ic(alpha2)).^2));

% minimize using fmincon
alpha2 = fmincon(obj2, alpha0, [], [], [], [], 0, [], '', opts);

orderA = alpha1(1);
orderB = alpha1(2);

ka = alpha1(2);
kb = alpha2(2);

% print stuff 
fprintf('The reaction order of reaction 1 is %0.0f \n',orderA)
fprintf('The reaction order of reaction 2 is %0.0f \n\n',orderB)

fprintf('The specific rate of reaction of reaction 1 is %0.2f \n',ka)
fprintf('The specific rate of reaction of reaction 2 is %0.2f \n\n',kb)

% plot stuff
figure(1)
hold on

scatter(t, Ca_im, 'b', 'o', 'LineWidth', 2);
plot(t, Ca_ic(alpha1), 'LineWidth', 2, 'color', 'c');

xlabel('Time (min)')
ylabel('Concentration (mol/dm^3)')
legend('Experimental','Model')
title('Concentration vs. Time for Reaction 1')

% plot stuff
hold off

figure(2)
hold on

scatter(t, Cb_im, 'r', 'o', 'LineWidth', 2);
plot(t, Cb_ic(alpha2), 'LineWidth', 2, 'color', 'm');

xlabel('Time (min)')
ylabel('Concentration (mol/dm^3)')
legend('Experimental','Model')
title('Concentration vs. Time for Reaction 2')

hold off

T  = [280 290 300 310 320 330 340 350];                 % give T
k1 = [0.306 0.408 0.510 0.660 0.816 1.030 1.21 1.36];   % given K1
k2 = [0.138 0.222 0.282 0.330 0.390 0.558 0.75 0.936];  % given K2

global R;
R  = 8.314;         % gas constant

reci_T = 1./T;      % reciprocal T
lnk1   = log(k1);   % ln(K1)
lnk2   = log(k2);   % ln(K2)

fit1 = polyfit(reci_T,lnk1,1);      % Fit line to lnK1 vs. 1/T 
fit2 = polyfit(reci_T,lnk2,1);      % Fit line to lnK2 vs. 1/T

yfit1 = polyval(fit1,reci_T);       % Fit line to lnK1 vs. 1/T 
yfit2 = polyval(fit2,reci_T);       % Fit line to lnK2 vs. 1/T 

global A1;
global A2;
    
global Ea1;
global Ea2;
    
Ea1 = -fit1(1)*R;   % -slope*R
Ea2 = -fit2(1)*R;   % -slope*R

A1 = exp(fit1(2));  % e^intercept
A2 = exp(fit2(2));  % e^intercept

% plot stuff
figure(3) 
hold on

scatter(reci_T, lnk1, 'b', 'o', 'LineWidth', 2)
scatter(reci_T, lnk2, 'r', 'o', 'LineWidth', 2) 

plot(reci_T, yfit1, 'LineWidth', 2, 'color', 'c')
plot(reci_T, yfit2, 'LineWidth', 2, 'color', 'm')

xlabel('1/Temperature (1/K)')
ylabel('ln(k)')
legend('Reaction 1','Reaction 2')
title('Arrhenius Plot for Reactions 1 and 2')

hold off

% print stuff
fprintf('The activation energy of reaction 1 is %0.0f Joules \n',Ea1)
fprintf('The activation energy of reaction 2 is %0.0f Joules \n\n',Ea2)

fprintf('The frequency factor of reaction 1 is %0.0f \n',A1)
fprintf('The frequency factor of reaction 2 is %0.0f \n\n',A2)


global Ca0 
global Cb0                       

global T0
global Tc
global Cpa
global Cpb

global delHa
global delHb   

Ca0 = 1;                      % mol/L
Cb0 = 2;                      % mol/L

T0  = 25+273;                 % initial temperature in K

Cpa = 25*4.184;               % Cp of A
Cpb = 50*4.184;               % Cp of B
   
delHa = 8000 *4.184;          % Delta H of A
delHb = 15000*4.184;          % Delta H of B    

Tc   = 5 + 273;               % Cooling water temperature in K
tab1 = funtime();

Tc   = 10 + 273;              % Cooling water temperature in K
tab2 = funtime();             

Tc   = 15 + 273;              % Cooling water temperature in K
tab3 = funtime();

disp(tab1)
disp(tab2)
disp(tab3)


function [tab] = funtime()
   
    % make a matirx (3 row 9 columns) to store all 
    % the values from fmincon
    table_vars = zeros(3,9);

    index = 1;

    global V

    for V = [100, 200, 500]

        A   = [];                     % Linear inequality constraints
        b   = [];                     % Linear inequality constraints

        Aeq = [];                     % Linear equality constraints
        beq = [];                     % Linear equality constraints

        lb  = [0   0.02 0.02 20 0   0   0   300];   % Lower bounds
        ub  = [inf inf  inf 100 100 inf inf 500];   % upper bounds

        opts = optimoptions(@fmincon, 'Display', 'off');
        [x, fval] = fmincon(@Profitfunction, [1 1 1 90 60 1 1 490], ...
                          A, b, Aeq, beq, lb, ub, @nonlincon, opts);
        fprintf('%0.10f \n\n',fval)

        % make table with each variable in a new column and a new row 
        % every iteration
        % round to make pretty display
        table_vars(index,1) = round(x(1),3);
        table_vars(index,2) = round(x(2),3);
        table_vars(index,3) = round(x(3),3);
        table_vars(index,4) = round(x(4),0);
        table_vars(index,5) = round(x(5),0);
        table_vars(index,6) = round(x(6),1);
        table_vars(index,7) = round(x(7),3);
        table_vars(index,8) = round(x(8),0);
        table_vars(index,9) = round(-fval,2);

        index = index+1;

    end

    % vector of all row names
    rowname = {'100 L'; '200 L'; '500 L'};

    % vector of all column names
    colname = {'Ca (mol/L)'; 'Cb (mol/L)'; 'Cc (mol/L)'; ...
               'Fa (L/min)'; 'Fb (L/min)'; 'mc (kg/hr)'; ...
               'H_fill (m)'; 'T (C)'; 'Profit ($/yr)'};

    % make table
    tab = table( ...
          table_vars(:,1), table_vars(:,2), table_vars(:,3), ...
          table_vars(:,4), table_vars(:,5), table_vars(:,6), ...
          table_vars(:,7), table_vars(:,8), table_vars(:,9), ...
          'VariableNames', colname, 'RowNames', rowname);

end


function [Profit] = Profitfunction(x)
   
    global Ca0                    
    global Cb0                       
    global V
    global Tc                     
    
    costA    = 50/1000;           % $/mol
    costB    = 20/1000;           % $/mol
    costC    = 200/1000;          % $/mol
    costR    = 6000*(V/1000)^0.8; % $
    costProc = 50/(100*1000);     % $/mol
  
    hrs_year = 8400;              % hrs of operation per year
    
    Ca = x(1);
    Cb = x(2);
    Cc = x(3);
    Fa = x(4);
    Fb = x(5);
    mc = x(6);
    
    if    Tc  == 278 
        costW = 1/100;            % $/kg
    
    elseif Tc == 283
        costW = 0.1/100;          % $/kg
    
    else 
        costW = 0.05/100;         % $/kg
    end
    
    Profit = -((hrs_year*60)*(Cc*(Fa+Fb)*costC - Ca0*Fa*costA ...
             - Cb0*Fb*costB - mc*costW/60 ...
             - (Ca + Cb + Cc)*(Fa + Fb)*costProc)  - costR);
end


function [C,Ceq] = nonlincon(x)
    
    global A1
    global A2
    
    global Ea1
    global Ea2
    
    global Ca0
    global Cb0
    global V
    
    global T0
    global Tc
    global Cpa
    global Cpb
    
    global R

    global delHa
    global delHb 
    
    
    Ca    = x(1);
    Cb    = x(2);
    Cc    = x(3);
    Fa    = x(4);
    Fb    = x(5);
    mc    = x(6);
    hfill = x(7);
    T     = x(8);
    
    % V = (pi*r^2*h) ==> pi*(D/2)^2*(2D) ==> (pi*D^3)/2
    % m 
    D = (2*(V/1000)/pi)^(1/3);  
    
    % 0.2*Vtank <= Vfill <= 0.9*Vtank
    % Vfill = pi*(D/2)^2*hfill ==> Vfill = pi*D^2/4*hfill 
    % ==> hfill = Vfill*4/(pi*D^2)
    % hfill = 4*Vfill/(pi*D^2) 
    % in m^3
    Vfill = (hfill*pi*D^2/4);  
    
    % mole balance of A
    A = Fa*Ca0 - (Fa+Fb)*Ca - A1*exp(-Ea1/(R*T))*Ca^2*(Vfill)*1000;
    
    % mole balance of B
    B = Fb*Cb0 - (Fa+Fb)*Cb + A1*exp(-Ea1/(R*T))*Ca^2*(Vfill)*1000 ...
        - A2*exp(-Ea2/(R*T))*Cb^2*(Vfill)*1000;
    
    % mole balance of C
    C = -(Fa+Fb)*Cc  + A2*exp(-Ea2/(R*T))*Cb^2*(Vfill)*1000;     
    
    pfill = (Vfill)/(V/1000);
    
    % A = pi*D*h 
    % multiply area by 0.75 since heat transfer only occurs from
    % 75% of the area
    UA = (0.75*pi*D*hfill)*(100*(mc)^0.8)*pfill*4.184;
    
    % energy balance 
    Energy = Fa*Cpa*Ca0*(T0-T) + Fb*Cpb*Cb0*(T0-T) - UA*(T-Tc) ...
             + delHa*A1*exp(-Ea1/(R*T))*Ca^2*Vfill*1000 ...
             + delHb*A2*exp(-Ea2/(R*T))*Cb^2*Vfill*1000;
    
    Ceq = [A; B; C; Energy];
    
    C   = [0.2*(V/1000)-(Vfill); (Vfill)-0.9*(V/1000)];
    
end
