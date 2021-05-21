clear all;
PreferenceParam = [10 0.9868 1.5] ;
StateParam=[0.00467 0.999 0.00062208 0.018]
FirmParam=[0.10 0.01 5];
global sigmav sigmax kappav vbar gamma delta psi theta a0 a1 mu
global k1 a b k0 d
global B kappa1
global  V nn uend W2 DN dt tt sigma x0
sigmav=StateParam(1); %volatility of the volatility process
sigmax=FirmParam(1); %idiosyncratic risk of the firm
kappav=StateParam(2); %Captures the speed of reversion of the volatility process 
vbar =StateParam(3); % Long run mean level of the process Vt
mu = StateParam(4); % Mean growth rate of the consumption process 
gamma= PreferenceParam(1); % Risk aversion parameter 
delta=PreferenceParam(2); % Subjective discount factor
psi=PreferenceParam(3); % Intertemporal elasticity of substitution 
theta=(1-gamma)/(1-(1/psi)); 
a0=FirmParam(2); % mean growth rate of the firm
a1=FirmParam(3); % systemic risk of the firm
T=10; % Maturity
d = 0; % d was needed to adjust for the boundary condition at v=0, but we do not need this anymore

% Solving the following two non-linear system of equations for kappa1 and B
% terms
kappa10=delta;
B0=((theta^2*kappa10^2*sigmav^2)^(-1))*(theta*(kappav+(1-kappa10))+sqrt((theta*(kappav+(1-kappa10)))^2-theta^4*kappa10^2*sigmav^2*(1-1/psi)^2));
fun = @BK;
x0 = [kappa10,B0];
y = fsolve(fun,x0);
B = y(2);
kappa1=y(1);

% rt is the risk-free interest rate given by rt = r0+r1Vt 
r1=(1-theta)*B*((kappa1-1)-kappa1*kappav)-(gamma^2+(1-theta)^2*kappa1^2*B^2*sigmav^2)/2;
r0=-(theta*log(delta)+(theta-1)*(-log(kappa1)-(kappa1-1)*B*vbar)-gamma*mu+...
    -(1-theta)*kappa1*B*kappav*vbar);

%the 
KK=kappav+(1-theta)*kappa1*B*sigmav^2;
TT=r1*kappav*vbar/(kappav+(1-theta)*kappa1*B*sigmav^2);
HH=sqrt(KK^2+2*r1*sigmav^2);

% k0+k1V_t is the drift of V_t with repect to the risk neutral measure
k1=-(kappav+(1-theta)*B*sigmav^2*kappa1);
k0=kappav*vbar-d*sigmav^2;

% Transformation constants a and b
%a=gamma/a1;
a=0;
b=k1*vbar/sigmav^2;
%b=0;

%the range for y
My=15/abs(b)
%% Coefficients of the pde 
f = 0;
m = 0;

%% Geometry
% Solve the problem on a square domain. The |squareg| function describes
% this geometry. Create a |model| object and include the geometry. Plot the
% geometry and view the edge labels.
numberOfPDE = 1;
model = createpde(numberOfPDE);
R1 = [3,4,-0.5,3,3, -0.5, My,My,-My,-My]';
gm = R1;
sf = 'R1';
ns = char('R1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);


%% Specify PDE Coefficients
% Define the function of coefficient 'a'
acoeffunction = @(region,state)r0+r1*(vbar*region.y+vbar) + (b/vbar)*(k0+k1*(vbar*region.y+vbar))+d*((kappav*vbar-sigmav^2)./region.y+k1)+...
                a*(a0+a1*gamma*region.y)-a^2*(a1^2*region.y+sigmax^2)-0.5*(b/vbar)^2*(sigmav^2*(vbar*region.y+vbar));
specifyCoefficients(model,'m',m,'d',1,'c',@ccoeffunction,'a',acoeffunction,'f',f);

%% Boundary Conditions
% Set zero Dirichlet boundary conditions on the left (edge 4), top (edge 1) and 
% bottom (edge 3) and Dirichlet boundary conditions on the right (edge 2).
applyBoundaryCondition(model,'dirichlet','Edge',[1  3 4],'u',0);

% Define the function of boundary condition on edge 2
myufunction = @(region,state)((2*HH*exp((HH+KK).*(state.time)/2)/(2*HH+(HH+KK)*...
               (exp(HH*(state.time))-1)))^(2*kappav*vbar/sigmav^2)).*... % At
               exp(-r0*state.time-(vbar*region.y+vbar).*(2*r1*(exp(HH*state.time)-1)./(2*HH+(HH+KK)*(exp(HH.*(state.time))-1)))).*... % Bt
              exp(a.*region.x+b.*region.y);
%myufunction = @(region,state)((2*HH*exp((HH+KK)*(state.time)/2)/(2*HH+(HH+KK)*...
               %(exp(HH*(state.time))-1)))^(2*kappav*vbar/sigmav^2))%.*... % At
               %exp(-(r0+r1*(vbar*region.y+vbar))*((exp(HH*state.time)-1)./(2*HH+(HH+KK)*(exp(HH*(state.time))-1))))... % Bt
              %.*exp(a*region.x+b*region.y);
applyBoundaryCondition(model,'dirichlet','Edge', 2 ,'u',myufunction,'Vectorized','on');

%% Generate Mesh
% Create and view a finite element mesh for the problem.
generateMesh(model,'Hmax',0.02);


%% Create Initial Conditions
u0=@(location) (1-(location.y ==-My)).*(1-(location.x == 3)).*(1-(location.x == -0.5)).*(1-(location.y ==My)).*exp(a*location.x+b*location.y).*(location.y).^d;
%(1-(location.x == 3)).*(1-(location.x == -0.5)).*(1-(location.y== 0.001)).*(1-(location.y ==1)).*exp(a*location.x+b*location.y).*(location.y).^d;
setInitialConditions(model,u0);

%% Define Solution Times
% Find the solution at n equally-spaced points in time from 0 to 10.
n = 1000;
tlist = linspace(0,10,n+1);

%% Calculate the Solution 
% Set the |SolverOptions.ReportStatistics| of |model| to |'on'|.
model.SolverOptions.ReportStatistics ='on';
result = solvepde(model,tlist);
u = result.NodalSolution;
%% Undo the transformation 
[p,e,t] = meshToPet(model.Mesh);
gbar=(exp(-a*p(1,:)-b*p(2,:)))';
g=gbar.*u;
%% Calculate yield for specific variance and specific firm value
v0=0;
N=n/10;
maturity=(10/N)*(1:N);
k2=find(p(2,:)>v0-0.01&p(2,:)<v0+0.01);
p1=p(1,k2);
[p1,I]=sort(p1);
g1=[];
x0=0.0665;
[ d, ix ] = min( abs( (p1')-x0 ) );
for i=1:N
   g1(i,:)=g(k2,i*10+1)';
   end
g1=g1(:,I);
gg=[];
for i=1:N
if abs(log(g1(i,ix)/g1(i,end)))< 0.0010
gg(1,i)=0;
else    gg(1,i)=-(log(g1(i,ix)/g1(i,end)))/maturity(1,i)
end;
end;
plot(maturity, gg(1,1:N));
%hold
%%%%
newx=[];
sigma=0.25;
DN=81; %DN-1 must be divisor of 2000 and be divisible by 4
rho=zeros(DN+1,526);
massrho=zeros(DN+1,1);
nrho=zeros(DN+1,526);
tt=linspace(0,1,DN);               
dt=1/(DN-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('trialsX.mat');%trialsX.mat contains a scenario for the firm value process between time 0 and 1 under the parameters set above.  If any of the parameters change, this file needs to be updated. 
load('trials'); %We use previously generated scenarios for volatility and consumption.  The data for these are in the data file "trial.mat"
trialnumbers=[9 10 23 24 31 32 33 34];
trn=1;
for trn=1:4
trn
W2=trial(1:2001, trialnumbers(trn*2))'; 
V=trial(1:2001, trialnumbers(trn*2-1));
X=trialX(1:2001,trn); 
%% Solve spde rho(t,x)
for nn=1:DN-1
m = 0;
t = linspace(tt(nn),tt(nn+1),101);
x = linspace(-0.5,10,526);
dx=(10-(-0.5))/526;
sol = pdepe(m,@rhofunction1pde,@rhofunction1ic,@rhofunction1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
uend=u(end,:)';           % value of rho at end point of each subintervel

if nn==1
     rho(1,:)=u(1,:);
     rho(nn+1,:)=u(end,:);
  else
      rho(nn+1,:)=u(end,:);
  end

end

newmean=log(X(2001)+0.5)+sigma^2;
newmean;
                      % at t=1, use the density below
rho(DN+1,:) = exp(-0.5*sigma^(-2)*(log(x+0.5)-newmean).^2).*((x+0.5+0.0001)*sigma*sqrt(2*pi)).^(-1);


newx=[newx, X(1500), X(2001)];
yieldimp=[];
for nn=1:DN+1
    massrho(nn,1)=sum(rho(nn,:)*dx);
   nrho(nn,:)=rho(nn,:)./massrho(nn,1);
   RF=[x', rho(nn,:)'];
   samplePoints = {x, 1:size(RF,2)};
   GI = griddedInterpolant(samplePoints,RF);
   if (nn < DN+1) vn=(V((nn-1)*2000/(DN-1)+1)-vbar)/vbar;
   else vn=(V(2001)-vbar)/vbar;
   end;
   kn=find(p(2,:)>vn-0.01&p(2,:)<vn+0.01);
pn=p(1,kn);
[pn,In]=sort(pn);
gn=[];
for i=1:N
   gn(i,:)=g(kn,i*10+1)';
end;
gn=gn(:,In);
TBprice=gn(:,end)',
BF=[pn', gn'];
[Q,m,l]=unique(BF(:,1));
BFn=BF(m,:);
xQ=unique(sort([x(1:176)';Q]));
queryPoints = {xQ,1:size(RF,2)};
   RFq = GI(queryPoints);
%interpolating B   
samplePointsb = {Q, 1:size(BFn,2)};
GIb = griddedInterpolant(samplePointsb,BFn);
queryPointsb = {xQ,1:size(BFn,2)};
BFq=GIb(queryPointsb);
dxQ=[0;diff(xQ)];
dxQ= repmat(dxQ,1,100);
RFq=repmat(RFq(:,2),1,100);
Price=[];
Price= RFq.*dxQ.*BFq(:,2:101);
Price=sum(Price)./sum(RFq.*dxQ);
yieldimp(nn,:)=-log(Price./TBprice)./maturity;
end;
h= figure;
for j=1:5
    plot(x,nrho(0.25*(DN-1)*(j-1)+1,:));
    axis([ -0.7 0.7 0 7.5]);
hold on
end;
plot(x,nrho(DN+1,:));
axis([ -0.7 0.7 0 12])
legend('t=0', 't=0.25', 't=0.50', 't=0.75', 't=1-', 't=1')
saveas(h,sprintf('rho%d.png',trn));
hold off

%Plotting the hazard rate of the default time
lambda=zeros(DN+1,1)
for i=1:DN
lambda(i)=0.5*(a1^2*V((i-1)*2000/(DN-1)+1)+sigmax^2)*(1/dx)*(nrho(i, 2)-nrho(i, 1))
end;
lambda(DN+1)= 0.5*(a1^2*V(2001)+sigmax^2)*(1/dt)*(nrho(DN+1, 2)-nrho(DN+1, 1))
h= figure;
plot(tt(1:DN-1),lambda(1:DN-1), '-*');hold on;
plot(tt,lambda(1:DN), '-o');hold on;
plot(1,lambda(DN+1), 'r*');hold on;
plot(1,lambda(DN+1), 'ro')
saveas(h,sprintf('lambda%d.fig',trn));
hold off

%plotting the yield spreads with imprefect information

h= figure;
for j=1:5
    plot(maturity,yieldimp(0.25*(DN-1)*(j-1)+1,:));
    axis([ 0 10 0 0.35]);
hold on
end;
plot(maturity,yieldimp(DN+1,:));
legend('t=0', 't=0.25', 't=0.50', 't=0.75', 't=1-', 't=1')
saveas(h,sprintf('yieldimp%d.fig',trn));
hold off




%% Calculate perfect yield, risk neutral yield at t=0.75 and compare to imperfect yield 
v3=(V(1500)-vbar)/vbar;
N=n/10;
maturity=(10/N)*(1:N);
k2=find(p(2,:)>v3-0.01&p(2,:)<v3+0.01);
p1=p(1,k2);
[p1,I]=sort(p1);
g1=[];
x3=X(1500);
[ d, ix ] = min( abs( (p1')-x3 ) );
for i=1:N
   g1(i,:)=g(k2,i*10+1)';
   end
g1=g1(:,I);
gg=[];
for i=1:N
if abs(log(g1(i,ix)/g1(i,end)))< 0.0010
gg(1,i)=0;
else    gg(1,i)=-(log(g1(i,ix)/g1(i,end)))/maturity(1,i)
end;
end;
plot(maturity, gg(1,1:N));
hold
plot(maturity,yieldimp(0.25*(DN-1)*3+1,:));
    axis([ 0 10 0 0.08]);
    saveas(h,sprintf('yieldimpperfect%d.fig',trn));
    hold off;
end