function u0 = rhofunction1ic(x)
  global uend nn sigma x0
  
  meanrho=log(x0+0.5)+sigma^2;
   if nn==1                     % at t=0, use x^2/(1+x^8) as initial value
      u0 = exp(-0.5*sigma^(-2)*(log(x+0.5)-meanrho)^2)*((x+0.5+0.0001)*sigma*sqrt(2*pi))^(-1);
   else
      u0=uend(round(x*50)+26);     % use values at end point as the initial value for next interval
   end
end