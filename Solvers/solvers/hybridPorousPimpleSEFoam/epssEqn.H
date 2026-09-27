fvScalarMatrix epssEqn
    (
       fvm::ddt(rhos,epss) 
       == -Ri
    );
epssEqn.solve();

forAll(epss,celli) 
    {

	if(epss[celli] >= 1) 
	{	
	   epss[celli] = 1; 
	}
        else if(epss[celli] <= 0.01) 
	{	
	   epss[celli] = 0; 
	}
    }
forAll(mesh.boundary(), patchi) 
    {
	forAll(epss.boundaryField()[patchi],facei)
	{
		if(epss.boundaryField()[patchi][facei] >= 1) 
		{epss.boundaryFieldRef()[patchi][facei] = 1;}
		else if (epss.boundaryField()[patchi][facei] <=0.01) 
		{epss.boundaryFieldRef()[patchi][facei] = 0;}
	}
     }


forAll(epss,celli) 
    {
	if(epss[celli] >= 0.01) 
	    {Solid[celli] = 1;}
        else
            {Solid[celli] = 0;}
    }
forAll(mesh.boundary(), patchi) 
    {
	forAll(epss.boundaryField()[patchi],facei)
	{
		if(epss.boundaryField()[patchi][facei] >= 0.01) 
		{Solid.boundaryFieldRef()[patchi][facei] = 1;}
		else
		{Solid.boundaryFieldRef()[patchi][facei] = 0;}
	}
     }
Solidf=fvc::interpolate(Solid);
forAll(Solidf,facei)
     {
	if(Solidf[facei]<1) 
	{Solidf[facei]=0;}
     }

eps=1-epss;

