function [J]= Jcalc(u,ud,dt)
% calculate the objective function J(eta1,eta2)

Udiff = u-ud;
Nt   = size(Udiff,2);
Ix  = zeros(1,1); % integration at constant time along x
for i = 1:Nt
    % abs is taken because (u-ud)^2
    I_loc  = Udiff(:,i)'*Udiff(:,i); % chebfun takes the integral
    Ix(i) = abs(I_loc);
end
J = sum(dt*(Ix(2:end)+Ix(1:end-1))/2);
end
