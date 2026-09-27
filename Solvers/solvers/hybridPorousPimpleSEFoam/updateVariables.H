k0=K*(eps*eps*eps)/(max((1-eps)*(1-eps),SMALL));

//Calculating Drag Coefficient
Drag= fvc::average(laminarTransport.nu()/k0)*Solid;


