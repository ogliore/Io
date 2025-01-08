clearvars
close all;

ntpoints=1e6;
Rot=logspace(1,3,1e2);
Ro=10e-6;

dynamic_viscosity=1;
rho_t=2800; % kg/m^3
nuL=dynamic_viscosity./rho_t; % m^2/s
sigmaL=0.35; % N/m
Pinf=0;

kvalue=1.29; % SO2

g_Io=1.796; % m/s2
rho_Io_crust=2792; % kg/m^3

eruption_speed=1;

dike_to_vent_area_ratio=1;
dike_speed=eruption_speed.*dike_to_vent_area_ratio;
max_depth=6.5e3; %6.5e3;

%tmax=max_depth./dike_speed;
%tmax=1e-1;


%depth=max_depth-tvec*dike_speed;

%P = rho_Io_crust .* g_Io .* depth;

% Pd=rho_Io_crust .* g_Io .* (max_depth-t*dike_speed);

% ((pvconst/P).^(1/gammavalue))^(1/3)


Rvo=0;

pvconst=(rho_Io_crust .* g_Io .* max_depth).*(Ro.^3).^kvalue;

p_G0=(rho_Io_crust .* g_Io .* max_depth)+2*sigmaL/Ro; %Pinf - P(1) + 2*sigmaL/Ro;

% A_SO2=
% Vapor_Pressure_SO2 = @(T) A_SO2-(B_SO2/(T+C_SO2)) .* T>A;


% p0 * V0^g = constant
% p0 * V0^g = constant

tmax=max_depth./dike_speed;


% T1 = @(R,Rv) -(3./(2.*R)).*Rv.^2;
% T2 = @(R,Rv) -(4.*nuL./R.^2).*Rv;
% T3 = @(R)    -(2*sigmaL)/(rho_t.*R.^2);
% T4 = @(R,P)   P./(rho_t.*R); %0; %(P-Pinf)./(rho_t.*R);
% T5 = @(R)     (1./R).*(p_G0./rho_t).*(Ro./R).^(3*k);
% Ra_ = @(R,Rv,P) T1(R,Rv) + T2(R,Rv) + T3(R) + T4(R,P) + T5(R);

odefun = @(t, y) [y(2); -(3./(2.*y(1))).*y(2).^2 +  -(4.*nuL./y(1).^2).*y(2) + -(2*sigmaL)/(rho_t.*y(1).^2) + (-(rho_Io_crust .* g_Io .* (max_depth-t*dike_speed)))./(rho_t.*y(1)) + (1./y(1)).*(((rho_Io_crust .* g_Io .* max_depth)+2*sigmaL/Ro)./rho_t).*(Ro./y(1)).^(3*kvalue)];
tspan = [0 tmax];
y0=[Ro, Rvo];
t=linspace(0,tmax,ntpoints);
[t, Y] = ode23s(odefun, t, y0);
Pv=rho_Io_crust .* g_Io .* (max_depth-t*dike_speed);

RtRo=Y(:,1)./Ro;

semilogx(Pv./1e6,RtRo,'LineWidth',3)
hold on
semilogx(Pv./1e6, ((pvconst./Pv).^(1/kvalue)).^(1/3)./Ro,'k--','LineWidth',3)
hold off
xlim([1e-3,40])
xlabel('Pressure (MPa)','FontSize',16)
ylabel('R/Ro','FontSize',16)
set(gca,'Linewidth',1.5,'FontSize',17.5)
saveas(gca,['RayleighPlesset_' num2str(Ro*1e6) 'um.png'])

save(['RayleighPlesset_' num2str(round(Ro*1e6)) 'um.mat']);

X=log(Pv./1e6);
Y=RtRo;


% R(1)=Ro;
% Rv(1)=Rvo;
% Ra(1)=Ra_(R(1),Rv(1),P(1));
% 
% deltat=tvec(2:end)-tvec(1:end-1);
% 
% 
% for ii=2:nt
% 
%     Rv(ii)=Rv(ii-1)+Ra(ii-1)*deltat(ii-1);
%     R(ii)=R(ii-1)+Rv(ii)*deltat(ii-1);
% 
%     Ra(ii)=Ra_(R(ii),Rv(ii),P(ii));
% 
% end
% 
% plot(tvec,R/Ro); %xlim([0,1e-7]); %ylim([0.1,10]);