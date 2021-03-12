% Finite Difference,SurfacePlot,And Partical diff
%% Finite Difference method
D = 3.5*10^-5; % cm^2/s
Cf = 50; % Micrometers
dr = 0.01; % Micrometers
dt= 0.1;   % s
r_vec = linspace(0,7,30); 
t_vec = 0:dt:30;
O_mat= zeros(length(r_vec),length(t_vec));
O_mat(1,:) = 0; % initial concentration left end boundary
O_mat(end,:)= Cf; % inital concentration right end

%%% Integration using Euler's and FDM 
for tdr =  1:length(t_vec)-1 % integrates to tdx+1
    for idr = 2:length(r_vec)-1
    O_mat(idr,tdr+1) = O_mat(idr,tdr) + D*dt/(dr^2)*(O_mat(idr+1,tdr) - 2*O_mat(idr,tdr) +O_mat(idr-1,tdr));
    end
end
%%% Plot
[rr,tt] = meshgrid(t_vec,r_vec);

mesh(rr,tt,O_mat)
xlabel('Time')
ylabel('Radial Distance um')
zlabel('O2 defused')
title('Diffusion of O2 Through Capillary')