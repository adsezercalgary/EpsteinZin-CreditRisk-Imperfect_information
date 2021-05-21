
function [c,f,s] = rhofunction1pde(x,t,u,DuDx)

global sigmax a0 a1 gamma W2 V nn dt tt DN

c = 2/(sigmax^2+a1^2*V(round((t-tt(1))*2000)+1));

f = DuDx;

s = (2*DuDx*(a0+a1*gamma*V(round((t-tt(1))*2000)+1)-a1*sqrt(V(round((t-tt(1))*2000)+1))...
     *(W2(nn*(2000/(DN-1))+1)-W2((nn-1)*(2000/(DN-1))+1))/dt))/(sigmax^2+a1^2*V(round((t-tt(1))*2000)+1)); %scaling of W2 should be 2000/(DN-1).  2001 
 %is the size of W1
end
