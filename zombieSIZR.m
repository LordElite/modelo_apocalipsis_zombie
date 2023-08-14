
[t,y] = ode45(@zombie,[0 30],[500; 0;1;0]);
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4))
grid
xlim([0 30])
ylim([0 500])
xlabel('Tiempo [días]')
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'color',[0 0 0])
ylabel('Población [miles de habitantes]')
hl=legend('Susceptibles','Infectados','Zombies','Removidos')
%hl=legend('Susceptibles','Zombies')
set(hl, 'TextColor','k', 'Color','w', 'EdgeColor','b')

function dydt = zombie(t,y)
Pi = 0;
Beta=95e-4;
Alpha = 5e-3;
Zeta = 1e-4;
Delta = 1e-4;
Ro= 1;


dydt = [Pi-Beta*y(1)*y(3)-Delta*y(1);
        Beta*y(1)*y(3)-Ro*y(2)-Delta*y(2);
        Ro*y(2)+Zeta*y(4)-Alpha*y(1)*y(3);
        Delta*y(1)+Delta*y(2)+Alpha*y(1)*y(3)-Zeta*y(4)];


end
