%% Wilson Methanol(1)/Water(2) @333.15K
% Methanol(1) & Water(2) @333.15K
clear all
P  = [19.953,39.223,42.984,48.852,52.784,56.652,60.614,63.998,67.924,70.229,72.832,84.562]; % Pressure in kPa
x1 = [0,0.1686,0.2167,0.3039,0.3681,0.4461,0.5282,0.6044,0.6804,0.7255,0.7776,1];
y1E= [0,0.5714,0.6268,0.6943,0.7345,0.7742,0.8085,0.8383,0.8733,0.8922,0.9141,1];
y2E= [1-y1E];
x2 = [1-x1];
T  = 335.15; % K
Zc1 = 0.224; % Methanol Critical Compressibility
Zc2 = 0.230; % Water Critical Compressibility
Tc1= 512.5;  % K
Tc2= 647.27; % K
Vl1= 114;    % cm3/gmol 
Vl2= 56;     % cm3/gmol
R  = 83.14;   % bar*cm3/molK
Tr1= T/Tc1;
Tr2= T/Tc2;
GammaE1 = zeros(1,length(P)); 
GammaE2 = zeros(1,length(P));
GammaR1 = zeros(1,length(P));
GammaR2 = zeros(1,length(P));
GammaMC1= zeros(1,length(P));
GammaMC2= zeros(1,length(P));
GEE     = zeros(1,length(P));
GEM     = zeros(1,length(P));
PbubbleE = zeros(1,length(P));
PdewE    = zeros(1,length(P));
PbubbleM = zeros(1,length(P));
PdewM    = zeros(1,length(P));
%% Solving For Saturated Pressure Using Antoine Eqn in kPa
Psat1 = (10^((7.89750)-((1474.08)/(229.13+(333.15-273.15)))))/7.5006;
Psat2 = (10^((8.01195)-((1698.785)/(231.04+(333.15-273.15)))))/7.5006;

%% Experimental Gamma Values
GammaE1(2:end-1) = ((y1E(2:end-1).*(P(2:end-1)))./((x1(2:end-1).*(Psat1))));
GammaE2(2:end-1) = ((y2E(2:end-1).*(P(2:end-1)))./((x2(2:end-1).*(Psat2))));

%% GE/RT Experimental Calculation 
GEE(2:end-1) = (x1(2:end-1).*log(GammaE1(2:end-1))) + (x2(2:end-1).*log(GammaE2(2:end-1)));

%% Solving For Saturated Liquid Volumes Using Racket Eqn
Vlsat1 = Vl1*Zc1^((1-Tr1)^(2/7)); % cm^3/gmol
Vlsat2 = Vl2*Zc2^((1-Tr2)^(2/7)); % cm^3/gmol
%% Optimization With OBJ Function
Fun  = @(h) (sum((((-x1(2:end-1).*(log(x1(2:end-1) + ((Vlsat2./Vlsat1).*exp(-h(1)./(R*T))).*x2(2:end-1)))) - (x2(2:end-1).*(log(x2(2:end-1) + ((Vlsat1./Vlsat2).*exp(-h(2)./(R*T))).*x1(2:end-1))))) - GEE(2:end-1))./GEE(2:end-1)).^2)./(length(P)-2); 
h0   = [5,5];
[Solu,OBJ,l] =fminsearch(Fun,h0);

%% Setting Wilson Constants
Alph1 = (Vlsat2/Vlsat1)*exp(-Solu(1)/(R*T));
Alph2 = (Vlsat1/Vlsat2)*exp(-Solu(2)/(R*T));
%% Calculating Gibbs Excess Model
GEM(2:end-1) = (-x1(2:end-1).*log(x1(2:end-1)+Alph1.*x2(2:end-1))) - (x2(2:end-1).*log(x2(2:end-1)+Alph2.*x1(2:end-1)));
%% Calculating Gamma Values
GammaR1 = -log(x1 + Alph1.*x2) + x2.*((Alph1./(x1+Alph1.*x2)) - (Alph2./(x2+Alph2*x1)));
GammaR2 = -log(x2 + Alph2.*x1) - x1.*((Alph1./(x1+Alph1.*x2)) - (Alph2./(x2+Alph2.*x1)));
GammaMC1 = exp(GammaR1);
GammaMC2 = exp(GammaR2);

%% Solving For Experimental Bubble Pressure
PbubbleE(2:end-1) = (x1(2:end-1).*GammaE1(2:end-1).*Psat1) + ((x2(2:end-1)).*GammaE2(2:end-1).*Psat2);
PbubbleE(1)       = Psat2;
PbubbleE(end)     = Psat1;

%% Solving For Vapor Mole Fraction 
y1E = (x1.*GammaE1.*Psat1)./(PbubbleE);
y2E = (1-y1E);

%% Solving For Experimental Dew Pressure
PdewE(2:end-1) = 1./((y1E(2:end-1)./Psat1)+(y2E(2:end-1)./Psat2));
PdewE(1)       = Psat2;
PdewE(end)     = Psat1;

%% Solving For Model Bubble Pressure
PbubbleM(2:end-1) = x1(2:end-1).*GammaMC1(2:end-1).*Psat1 + x2(2:end-1).*GammaMC2(2:end-1).*Psat2;
PbubbleM(1)       = Psat2;
PbubbleM(end)     = Psat1;

%% Solving For Model Vapor Mole fractions
y1M = (x1.*GammaMC1.*Psat1)./(PbubbleM);
y2M = 1-y1M;

%% Solving For Model Dew Pressure
PdewM(2:end-1) = 1./((y1M(2:end-1)./Psat1)+(y2M(2:end-1)./Psat2));
PdewM(1)       = Psat2;
PdewM(end)     = Psat1;

%% Plotting
subplot(2,1,1)
scatter(x1,GEE)
hold on
plot(x1,GEM)
title('Methanol(1)/Water(2)')
xlabel('Methanol')
ylabel('GE/RT')
legend('Experimental(Raoults)','Model')
hold off

subplot(2,1,2)
hold on
scatter(x1,PbubbleE)
scatter(x1,PdewE)
plot(x1,PbubbleM)
plot(x1,PdewM)
xlabel('Methanol')
ylabel('Pressure(kPa)')
legend('Experimental Bubble', 'Experimental Dew','Model Bubble','Model Dew')
hold off

fprintf('Here we see the model fits the data well')
