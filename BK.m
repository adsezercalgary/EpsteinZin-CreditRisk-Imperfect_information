function y = BK(x)
global sigmav sigmax kappav vbar gamma delta psi theta a0 a1 mu

y(1) = theta*log(delta)-theta*log(x(1))+kappav*theta*x(1)*x(2)*vbar + theta*(1-x(1))*x(2)*vbar+ mu*(1-gamma);
y(2) = (theta^2*x(1)^2*sigmav^2)*x(2)^2-(theta*(kappav+(1-x(1))))*x(2)+0.5*theta^2*(1-1/psi)^2


