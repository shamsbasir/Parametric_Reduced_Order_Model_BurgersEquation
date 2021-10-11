function [u]=SolveState(t,u0,nu,eta)
dom = [-1 1];
x = chebfun('x',dom);
opts     = pdeset('Eps', 1e-4, 'Ylim', [-40,40]);
pdefun   = @(t,x,u) -u.*diff(u) + nu.*diff(u,2) +f(t,x,eta);
bc.left  = {'dirichlet', 0};
bc.right = {'dirichlet', 0};
u        = pde15s(pdefun,t,u0,bc,opts);
end
