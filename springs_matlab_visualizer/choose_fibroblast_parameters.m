
D0 = 1 ;
Dsat = 3.5173 ;

strain_min = 0.36 ;
strain_max = 0.43 ;
strain_delta = strain_max-strain_min ;

strain = 0.5*(strain_min+strain_max) ;
f = 0.2 ;
% ser = 0.5*(pi*f*(strain_delta))^2
ser = (1^2)*1*f*((strain_max^2)-(strain_min^2))

D = @(s) D0*exp(-10.4350*s) + Dsat./(1+exp(-13.0130*(s-0.3709))) ;

Dstar = D(strain)
Dmin = D( fmincon(D,0.5,[],[],[],[],0.0,1.0) )

R0 = Dmin ;
beta = 1/(2*ser)

Rstar = Dstar ;

Rsat = ((Rstar-R0)/(1-exp(-beta*ser))) + R0

R = @(e) R0 + (Rsat-R0)*(1-exp(-beta*ser)) ;
R( ser )
