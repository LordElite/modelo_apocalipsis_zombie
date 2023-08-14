

[t,y] = ode45(@zombie,[0 5],[500; 1;0]);
plot(t,y(:,1),t,y(:,2),t,y(:,3))
grid
xlim([0 5])
ylim([0 600])
xlabel('Tiempo [días]')
ylabel('Población [miles de habitantes]')
%set(gca,'Xcolor','w');
%set(gca,'Ycolor','w');
%set(gca,'color',[0 0 0])
hl=legend('Susceptibles','Zombies','Removidos')
set(hl, 'TextColor','k', 'Color','w', 'EdgeColor','b')

function dydt = zombie(t,y)
Pi = 0;
Beta = 95e-4;
%Beta=190e-4;
Alpha = 5e-3;
Zeta = 1e-4;
Delta = 1e-4;


dydt = [Pi-Beta*y(1)*y(2)-Delta*y(1);
    Beta*y(1)*y(2)+Zeta*y(3)-Alpha*y(1)*y(2);
    Delta*y(1)+Alpha*y(1)*y(2)-Zeta*y(3)];

end
