Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid velocity field 'U'\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid fraction field epss\n" << endl;
volScalarField epss
(
    IOobject
    (
        "epss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1-epss
);

surfaceScalarField nEpsf
(
   IOobject
   (
      "nEpsf",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("nEpsf",dimArea, 0.0)
);

// Defenition of Solid Indicator
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceVectorField interpolatedU = fvc::interpolate(U);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
singlePhaseTransportModel laminarTransport(U, phi);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar rhos("rhos",transportProperties.lookup("rhos"));
dimensionedScalar rho("rho",transportProperties.lookup("rho"));

//
Info<< "Reading SoildProperties\n" << endl;
IOdictionary SoildProperties
(
    IOobject
    (
        "SoildProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar thickness("thickness",SoildProperties.lookup("thickness"));
dimensionedScalar Dp("Dp",SoildProperties.lookup("Dp"));
dimensionedScalar dp("dp",SoildProperties.lookup("dp"));
dimensionedScalar kB("kB",SoildProperties.lookup("kB"));
dimensionedScalar Tw("Tw",SoildProperties.lookup("Tw"));
dimensionedScalar h("h",SoildProperties.lookup("h"));
dimensionedScalar AH("AH",SoildProperties.lookup("AH"));
dimensionedScalar F("F",SoildProperties.lookup("F"));
dimensionedScalar z("z",SoildProperties.lookup("z"));
dimensionedScalar Epsilon0("Epsilon0",SoildProperties.lookup("Epsilon0"));
dimensionedScalar Epsilonr("Epsilonr",SoildProperties.lookup("Epsilonr"));
dimensionedScalar sigma0("sigma0",SoildProperties.lookup("sigma0"));
dimensionedScalar R("R",SoildProperties.lookup("R"));
dimensionedScalar c("c",SoildProperties.lookup("c"));
dimensionedScalar ker("ker",SoildProperties.lookup("ker"));

dimensionedScalar sp("sp",3.1415926*dp*dp/4);
dimensionedScalar kappa("kappa",sqrt(2*F*F*c*z*z/(Epsilon0*Epsilonr*R*Tw)));
dimensionedScalar s0("s0",(z*F*sigma0)/(Epsilon0*Epsilonr*R*kappa*Tw));
dimensionedScalar yInfinity0("yInfinity0",2*asinh(s0/2));
dimensionedScalar yInfinityh("yInfinityh",4*atanh(tanh(yInfinity0/4)*exp(-kappa*h)));
dimensionedScalar yInfinitym("yInfinitym",4*atanh(tanh(yInfinity0/4)*exp(-kappa*h/2)));
dimensionedScalar ym("ym",asinh(2*sinh(yInfinitym)+(4/(kappa*h))*sinh(yInfinityh/2)));
dimensionedScalar FA("FA",(AH*sp/(6*3.1415926))*((pow(h,-3))-(2*pow((h+thickness),-3))+pow((h+2*thickness),-3)));
dimensionedScalar FB("FB",2*c*R*Tw*sp*(cosh(ym)-1));

dimensionedScalar abc("abc", dimensionSet(1,-1,-2,0,0,0,0),1);

volScalarField tauc
(
    IOobject
    (
        "tauc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        
    ),
    Solid*((FA+FB)/sp+kB*Tw/(6*3.1415926*h*h*h))
);

volScalarField tauf
(
    IOobject
    (
        "tauf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Solid*6*3.1415926*Dp*laminarTransport.nu()*rho*fvc::average(mag(interpolatedU))/sp
);

volScalarField def("def", max((tauf-tauc),0*abc));

volScalarField er
(
    IOobject
    (
        "er",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    ker*def
);

Info << "Reading reaction rate if present\n" << endl;
volScalarField Ri
(
    IOobject
    (
        "Ri",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ri",dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);

dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);



//
Info<< "Calculating averaged mass flux\n" << endl;
surfaceScalarField phiByEps
(
	IOobject
	(
        "phiByEps",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*linearInterpolate(1/(eps+SMALL))
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());

#include "createPorousMediaFields.H"

