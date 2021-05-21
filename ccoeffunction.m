function cmatrix = ccoeffunction(region,state)
global sigmav sigmax kappav a0 a1 kappa1 B vbar psi delta gamma theta k1 a d k0

%k0=kappav*vbar-d*sigmav^2
n1 = 4;
nr = numel(region.x);
cmatrix = zeros(n1,nr);
cmatrix(1,:) = (a1^2*(vbar*region.y+vbar)+sigmax^2)/2;
cmatrix(2,:) =(a0-a1*gamma*vbar-a*sigmax^2)*region.y-0.5*a1*(gamma-a*a1)*vbar*(region.y).^2-(k0/vbar-0.5*sigmav^2/vbar)*region.x;
%cmatrix(2,:) = (a0-gamma*sigmax^2/a1)*region.y+ 0.5*(k1)*(region.y).^2-(k0-0.5*sigmav^2)*region.x;
cmatrix(3,:) = -cmatrix(2,:);
cmatrix(4,:) = sigmav^2/vbar*(region.y+1)/2;
end
