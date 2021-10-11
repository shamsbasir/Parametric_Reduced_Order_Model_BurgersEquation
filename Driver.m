% Shamsulhaq Basir
clear ;  close  ; clc
LW = 'linewidth'; FS = 'fontsize'; IN = 'interpret'; LT = 'latex';
%% ------------- Space-Time Domains -----------

Tf  = 2.0; dt = 0.01;
t   = 0.0:dt:Tf;
nu  = 0.05;
dom = [-1 1];
%% ---------------- Initial state ------------
x   = chebfun(@(x) x, dom);
% Initial condition.
u0  = sin(pi*x);

%% --------------- POD Computation ---------
n = 5;
xi1      = linspace(-1,1,n); xi2= linspace(-1,1,n);
[Xi1,Xi2] = meshgrid(xi1,xi2);
A =[];
for j=1:n
    for i=1:n
        
        eta(1,1) = Xi1(j,i);
        eta(1,2) = Xi2(j,i);
        u        = SolveState(t,u0,nu,eta);
        A        = [A u];
        
    end
end
eta_d(1,1) = -0.75;
eta_d(1,2) =  0.25;
u_d = SolveState(t,u0,nu,eta_d);
[U,S,V] = svd(A);

sigma = diag(S);
sigma = sigma(1:25);
figure
semilogy(sigma,'-o',LW,1.2)
grid on;
set(gca,FS,12)
ylabel("Singular value, $\sigma_k$",IN,LT);
xlabel("r");



figure()
plot(cumsum(sigma)/sum(sigma)*100,'-o',LW,1.2)
grid; set(gca,FS,12)
ylabel("Cumulative Energy");
xlabel("r");

for i = 1:numel(sigma)
    percentage = 100 - cumsum(sigma(i))/sum(sigma)*100;
    if(percentage >= 99.99)
        r = i;
        break;
    end
end

%-- 3 POD modes
figure()
plot(U(:,1));
hold on
plot(U(:,2));
plot(U(:,3));
xlabel('X');
ylabel('$\Phi(x)$');
legend('mode 1','mode 2','mode 3');

save data


%% -------------------------------------------------------------------

load data
%% (4) Response Surface
Tf = 2.0; dt = 0.01;
t = 0.0:dt:Tf;
snap = numel(t);
nu = 0.05;
dom = [-1 1]; x = chebfun('x',dom);
u0 = sin(pi*x);

Ubase   = U(:,1:r);                             % Basis
Ubase_p = diff(Ubase);                          % Derivative of Basis
Mik     = Ubase'*Ubase;                         % Mass matrix = I
y0      = (u0'*Ubase);                          % projection coefficients, t = 0
y0      = y0';
Dik     = -Ubase_p'*Ubase_p;                    % D matrix


dimEta = 9;
xi1=linspace(-1,1,dimEta);
xi2=linspace(-1,1,dimEta);
[Xi1,Xi2]= meshgrid(xi1,xi2);
J = zeros(1,1);

for j=1:dimEta
    for i=1:dimEta
        
        eta(1,1) = Xi1(i,j);
        eta(1,2) = Xi2(i,j);
        
        dydt = @(t,y) nu*Dik*y-Ubase'*((Ubase*y).*(Ubase_p*y))+Ubase'*f(t,x,eta);
        [t,y]=ode23(@(t,y) dydt(t,y), t, y0);
        u_rom = Ubase*y';                    % reduced order model solution
        J(i,j) =Jcalc(u_rom, u_d,dt);
    end
end
% find the minimum of J
[nx,ny] = find(J == min(min(J)));

%% Response surface
figure
surf(Xi1,Xi2,J)
et1 = Xi1(nx,ny);
et2 =Xi2(nx,ny);
hold on
plot3(et1,et2,J(nx,ny),'^','color','r');
xlabel('$\zeta_{1}$',FS,12);
ylabel('$\zeta_{2}$',FS,12);
zlabel('J($\zeta_{1}$ ,$\zeta_{1}$)',FS,12);
set(gca,FS,12);
%% optimum point
disp("Optimum is : ");
fprintf("(%2.3f, %2.3f) \n", Xi1(nx,ny),Xi2(nx,ny));

