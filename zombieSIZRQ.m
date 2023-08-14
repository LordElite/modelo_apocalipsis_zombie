
[t,y] = ode45(@zombie,[0 30],[500; 0;1;0;0]);


plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4),t,y(:,5))
grid
xlim([0 30])
ylim([0 500])
%title(sprintf('\kappa= %f, \sigma= %f, \zeta= %f, \gamma= %f', Kappa, Sigma, Zita, Gamma))
xlabel('Tiempo [días]')
ylabel('Población [miles de habitantes]')
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'color',[0 0 0])
hl=legend('Susceptibles','Infectados','Zombies','Removidos','En cuarentena')

set(hl, 'TextColor','k', 'Color','w', 'EdgeColor','b')

function dydt = zombie(t,y)

N=500;
Pi=0;
Beta = 95e-4;
Alpha = 5e-3;
Zeta = 1e-4;
Delta = 1e-4;
Ro= 1;
Kappa = 1e-2;
Sigma = 5e-3;
Gamma = 0.5;
dydt = [Pi-Beta*y(1)*y(3)-Delta*y(1);
        Beta*y(1)*y(3)-Ro*y(2)-Delta*y(2)-Kappa*y(2);
        Ro*y(2)+Zeta*y(4)-Alpha*y(1)*y(3)-Sigma*y(3);
        Delta*y(1)+Delta*y(2)+Alpha*y(1)*y(3)-Zeta*y(4)+Gamma*y(5);
        Kappa*y(2)+Sigma*y(3)-Gamma*y(5)];


end