close all
clear
clc
%% SAW Propagation in Grating Structure
c1 = 2700; % Wave speed in gap
c2 = 340; % Wave speed in grating
p = 0.1; % Length of periodicity
w = 0.09; % Length of grating
R1 = 1+0.01*pi; % Acoustic impedance in grating
R2 = 1;  % Acoustic impedance in gap
Ce = 0; % Equivalent capacitance
N = 100; % Number of grating cell

f0 = c1*c2/(2*w*c2+2*(p-w)*c1); % Center frequency
dt = 0.03;
f = (1-dt)*f0:1:(1+dt)*f0; % Frequency range

beta1 = 2*pi.*f/c1;
beta2 = 2*pi.*f/c2;
theta1 = beta1.*w;
theta2 = beta2.*(p-w);
ome = 2*pi.*f; % Angular frequency

% Exact solution
% f11 = cos(theta1./2).*cos(theta2./2)-(R2/R1).*sin(theta1./2).*sin(theta2/2)-ome.*Ce.*R2.*cos(theta1./2).*sin(theta2./2);
% f22 = cos(theta1./2).*cos(theta2./2)-(R1/R2).*sin(theta1./2).*sin(theta2/2)-ome.*Ce.*R1.*sin(theta1./2).*cos(theta2./2);
% f12 = 1i.*R2.*(cos(theta1./2).*sin(theta2./2)+(R1/R2).*sin(theta1./2).*cos(theta2/2)-ome.*Ce.*R1.*sin(theta1./2).*sin(theta2./2));
% f21 = 1i.*(1/R2).*(cos(theta1./2).*sin(theta2./2)+(R2/R1).*sin(theta1./2).*cos(theta2/2)+ome.*Ce.*R2.*cos(theta1./2).*cos(theta2./2));
% thetae = acos(f11.*f22+f12.*f21);
% Re = sqrt(f11.*f12./f22./f21);

% Approximate solution Assuming R1~=R2 & ome*Ce*R2<<1
Delta = R1/R2-1;
Phi = ome.*(Ce*R2);
T = sign(Delta).*sqrt(Delta.^2+Phi.^2);
ita = atan(Phi./Delta);
thetae = acos((1+T.^2/4).*cos(theta1+theta2+Phi)-(T.^2/4).*cos(theta1-theta2-2.*ita));
Re = R2.*sqrt((sin(theta1+theta2+Phi)+T.*sin(theta1+Phi./2-ita))./(sin(theta1+theta2+Phi)-T.*sin(theta1+Phi./2-ita)));

ref = (theta1+theta2)./pi-1; % Normalized frequency

Gam0 = (Re-R2)./(Re+R2);
Gam = Gam0.*(1-exp(-2*1i*N.*thetae))./(1-(Gam0.^2).*exp(-2*1i*N.*thetae)) ;

%% Plot figures
figure(1)

subplot(1,2,1)
plot(ref,10.*log10(abs(Gam0)))
xlim([-dt dt])
xlabel('Relative frequency');
ylabel('Reflection coefficient in dB');
axis square
title('$\Gamma_{0}$','Interpreter','latex')
hold on

subplot(1,2,2)
plot(ref,180.*angle(Gam0)./pi)
xlim([-dt dt])
ylim([0 180])
xlabel('Relative frequency');
ylabel('Reflection coefficient angle');
axis square
title('$\Gamma_{0}$','Interpreter','latex')
hold on

angGam=180.*angle(Gam)./pi;
for i = 1:length(angGam)
    if (angGam(i)<-90)
        angGam(i)=angGam(i)+360;
    end
end

figure(2)

subplot(1,2,1)
plot(ref,10.*log10(abs(Gam)))
xlim([-dt dt])
xlabel('Relative frequency');
ylabel('Reflection coefficient in dB');
axis square
title('$\Gamma$','Interpreter','latex')
hold on

subplot(1,2,2)
plot(ref,angGam)
xlim([-dt dt])
ylim([-90 270])
xlabel('Relative frequency');
ylabel('Reflection coefficient angle');
axis square
title('$\Gamma$','Interpreter','latex')
hold on