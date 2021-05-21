clear all
StateParam=[0.00467 0.999 0.00062208 0.018];
sigmav=StateParam(1); %volatility of the volatility process
kappav=StateParam(2); %Captures the speed of reversion of the volatility process 
vbar =StateParam(3); % Long run mean level of the process Vt
mu = StateParam(4); % Mean growth rate of the consumption process 
trial=[];
for m=1:20
N=201;
t1=linspace(0,1,(N-1)*10+1);
dt1=t1(2)-t1(1);
m
V11=zeros((N-1)*10+1,1);                     % Volatility process Vt               
 V11(1)=vbar;
 dwtv=normrnd(0,sqrt(dt1),1,(N-1)*10);
 for i=2:(N-1)*10+1
     V11(i)=V11(i-1)+ kappav*(vbar-V11(i-1))*dt1+sigmav*sqrt(V11(i-1))*dwtv(i-1);
 end;
h=figure;
subplot(1,2,1);
plot(t1, V11);
ylim([0.0003 0.0008]);
xlabel('time'); 
ylabel('$V_t$', 'Interpreter', 'latex'); 
%make here your second plot
 W1=zeros((N-1)*10+1,1); % Generate a Brownian Wc, mu = 0, sigma^2=dt1
 lnC= zeros((N-1)*10+1,1);
dwt1=normrnd(0,sqrt(dt1),1,(N-1)*100);

for i=2:(N-1)*10+1
    W1(i)=W1(i-1)+dwt1(i-1);
    lnC(i)=lnC(i-1)+mu*dt1+ sqrt(V11(i-1))*dwt1(i-1);
end;
subplot(1,2,2);    
plot(t1, lnC);
ylim([-0.05 0.10])
xlabel('time') 
ylabel('$\ln C_t-\ln C_0$', 'Interpreter', 'latex') 
trial=[trial, V11,W1];

saveas(h,sprintf('scenario%d.png',m)); % will create FIG1, FIG2,...

end