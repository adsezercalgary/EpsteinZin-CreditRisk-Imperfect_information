clear all
StateParam=[0.00467 0.999 0.00062208 0.018]
FirmParam=[0.10 0.01 5];
a0=FirmParam(2); % mean growth rate of the firm
a1=FirmParam(3);
mu = StateParam(4);
sigmax=FirmParam(1);
Nx=201
t1=linspace(0,1,(Nx-1)*10+1);
dt1=t1(2)-t1(1);
Wx=zeros((Nx-1)*10+1,1); % Generate a Brownian Wx, mu = 0, sigma^2=dt1
dwtx=normrnd(0,sqrt(dt1),1,(Nx-1)*100);
trialX=[];
x0=0.0665;
load('trials');
trialnumbers=[9 10 23 24 31 32 33 34];
trn=1;
for trn=1:4
trn
X=zeros((Nx-1)*10+1,1)
W2=trial(1:2001, trialnumbers(trn*2))';
V=trial(1:2001, trialnumbers(trn*2-1));
X(1)=x0
for i=2:(Nx-1)*10+1
    Wx(i)=Wx(i-1)+dwtx(i-1);
    X(i)=X(i-1)+(a0+a1*mu)*dt1+ a1*sqrt(V(i-1))*(W2(i)-W2(i-1))+sigmax*dwtx(i-1);
end;
h=figure;
plot(t1, X);
xlabel('time'); 
ylabel('$X_t$', 'Interpreter', 'latex'); 
saveas(h,sprintf('Xscenario%d.png',trn));
trialX=[trialX,X]
end