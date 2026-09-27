interpolatedU = fvc::interpolate(U); 

tauc=Solid*((FA+FB)/sp+kB*Tw/(6*3.1415926*h*h*h));

tauf=Solid*6*3.1415926*Dp*laminarTransport.nu()*rho*fvc::average(mag(interpolatedU))/sp;

def=max((tauf-tauc),0*abc);

er=ker*def;

