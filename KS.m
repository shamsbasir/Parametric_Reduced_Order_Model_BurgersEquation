clc; close all; clear
% Problem Parameter
Tf  = 2.0; dt = 0.01; 
t   = 0.0:dt:Tf;
nu  = 0.05;
dom = [-1 1];
x = chebfun('x',dom);

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-4, 'Ylim', [-40,40]);
u0 = sin(pi*x);
xi1=linspace(-1,1,5); xi2=linspace(-1,1,5);
[Xi1,Xi2]= meshgrid(xi1,xi2);
A =[];
for j=1:5
    for i=1:5
 
        pdefun = @(t,x,u) -u.*diff(u) + nu.*diff(u,2) + f(t,x,Xi1(j,1),Xi2(j,i));
        bc.left = {'dirichlet', 0};
        bc.right = {'dirichlet', 0};
        u = pde15s(pdefun,t,u0,bc,opts);
        
        A = [A u];

    end
end
save('AData','A');

% solving for desired state
d_A =[];
pdefun = @(t,x,u) -u.*diff(u) + nu.*diff(u,2) - 0.75*sin(pi*x-t) + 0.25*sin(2*pi*x-3*t);
bc.left = {'dirichlet', 0}; bc.right = {'dirichlet', 0};
u = pde15s(pdefun,t,u0,bc,opts);
d_A = [d_A u];

save('desired_state','d_A');
