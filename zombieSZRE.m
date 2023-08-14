

N = 500;
tmax =10;
dt = 0.05;
Pi = 0;
Beta = 195e-4;
%Beta=190e-4;
Alpha = 5e-3;
Zeta = 1e-4;
Delta = 1e-4;
t_oleada =[3 5 7]; 
k = 0.25; 
[S,R,Z,t] = erad(N,Pi,Alpha,Beta,Zeta,Delta,k,t_oleada,tmax,dt);

clf;close all;
figure
plot(t,S,t,Z,t,R);
grid
xlim([0 10])
ylim([0 600])
axis tight
xlabel('Tiempo [días]')
ylabel('Población [miles de habitantes]')
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'color',[0 0 0])
hl=legend('Susceptibles','Zombies','Removidos')
set(hl, 'TextColor','k', 'Color','w', 'EdgeColor','b')

function [S,R,Z,t] = erad(Npop,birthrate,alpha,beta,zeta,delta,killRate,tRaid,tmax,dt)


N = round(1+tmax/dt); % pasos de tiempo
Y = zeros(3,N);
Y(:,1) = [Npop;0;0]; % SIZR

t = (0:1:N-1)*dt; % tiempo
modelFun = @(Y,A,F) A*Y + F;


for ii=1:N-1
    A = getA(zeta,delta);
    SZ = Y(1,ii)*Y(2,ii);
    F = [birthrate-beta;beta-alpha;alpha].*SZ;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
    
    if min(abs(t(ii)-tRaid))<= 0.1
        Y(2,ii+1) =  max(0,Y(2,ii) -killRate*Y(2,ii));
        
        [ ~,ind] = (min(abs(t(ii)-tRaid)));
        tRaid(ind)=[];
    end
    
end

t = (0:1:N-1)*dt; % tiempo
S = Y(1,1:N);
Z = Y(2,1:N);
R = Y(3,1:N);


    function [A] = getA(zeta,delta)
        S0 = [-delta,0,0];
        Z0 = [0,0,zeta];
        R0 = [delta,0,-zeta,];
        A = [S0;Z0;R0];
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        
        % Runge-Kutta de cuarto orden
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % resultado
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end



end



