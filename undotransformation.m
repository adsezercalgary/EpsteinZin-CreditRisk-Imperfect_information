%% Undo transformation
[p,e,t] = meshToPet(model.Mesh);
gbar=(exp(-a*p(1,:)-b*p(2,:))./p(2,:))';
g=gbar.*u;
figure
gmax = max(max(g));
gmin = min(min(g));
for i = 1:n
    pdeplot(model,'XYData',g(:,i),'ZData',g(:,i),'ZStyle','continuous',...
                  'Mesh','off','XYGrid','on','ColorBar','on');
    axis([-0.5 3 0 1 gmin gmax]); 
    caxis([gmin gmax]);
    xlabel x
    ylabel y
    zlabel u
    M(i) = getframe;
end