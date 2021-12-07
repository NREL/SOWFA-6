// Read ABLProperties with input data
IOdictionary ABLProperties
(
    IOobject
    (
        "ABLProperties",
        runTime.time().constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


//- Automatically update fixesValue for timeVaryingMappedInletOutlet BC?
bool updateFixesValue(ABLProperties.lookupOrDefault<bool>("updateFixesValue",false));


  // This switch dictates whether or not the pressure reference cell is set explicitly
  // in the p_rgh solve.  If true, pressure is constrained at the pressure reference
  // cell by manipulating the matrix row for that cell to remove the null space.  This
  // assures that the pressure level is constrained, but it also means the continuity
  // equation is no longer solved at that cell, so the divergence error can be significant
  // enough there to cause localized oscillations.  Iterative solvers can still converge
  // even with a null space, but more iterations may be needed, so this switch can be
  // set to false.
  bool activatePressureRefCell(pimple.dict().lookupOrDefault<bool>("activatePressureRefCell", true));
  if (activatePressureRefCell)
  {
       Info << "Pressure reference cell matrix row modification enabled" << endl;
  }
  else
  {
       Info << "Pressure reference cell matrix row modification NOT enabled" << endl;
  }



// Create and read the potential temperature field [K].
Info << "Creating and reading potential temperature field, T" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the modified pressure field [m^2/s^2].
Info << "Creating and reading the density-normalized modified pressure field, p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the velocity field [m/s].
Info << "Creating and reading velocity field, U" << endl;
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


// Create and calculate the velocity flux [m^3/s].
Info << "Creating and calculating velocity flux field, phi" << endl;
#include "createPhi.H"


// Set up the laminar (molecular) transport model.
Info << "Creating and reading laminar transport properties" << endl;
singlePhaseTransportModel laminarTransport(U, phi);


// Read in reference temperature [K] and laminar and turblent Prandtl numbers.
// Note that turbulent Prandtl number may be overriden by a stability-
// dependent value based on the choice of SGS turbulence model.
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);
dimensionedScalar Pr("Pr", dimless, laminarTransport);
dimensionedScalar Prt("Prt", dimless, laminarTransport);


// Create the turbulence model.
Info << "Creating and reading turbulent transport properties" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


// Create Coriolis force object [m/s^2].
Info << "Creating Coriolis force" << endl;
CoriolisForce Coriolis(U);


// Create external source term objects [m/s^2 and K/s].
Info << "Creating and reading volumetric source terms for momentum and temperature" << endl;
drivingForce<vector> momentumSourceTerm(ABLProperties, "momentum",U);
drivingForce<scalar> temperatureSourceTerm(ABLProperties, "temperature",T);


// Create list of sponge layer object
  Info << "Creating sponge layer(s)" << endl;
  spongeLayer spongeLayers("spongeList",U);


// Create and read the kinematic turbulent thermal conductivity [m2/s].
Info << "Creating and reading turbulent thermal conductivity field, kappat" << endl;
volScalarField kappat
(
    IOobject
    (
        "kappat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the wall shear stress field [m^2/s^2].
Info << "Creating and reading the wall shear stress field, Rwall" << endl;
volSymmTensorField Rwall
(
    IOobject
    (
        "Rwall",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the wall temperature flux field [K-m/s].
Info << "Creating and reading the wall temperature flux field, qwall" << endl;
volVectorField qwall
(
    IOobject
    (
        "qwall",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Create and read the density-normalied pressure field [m^2/s^2].
Info << "Creating and reading the density-normalized pressure field, p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh
);


// Deal with reference pressure level.
Info << "Setting reference pressure level" << endl;
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);


// Create buoyancy model, which requires a reference height, which is the height of
// the pressure reference cell.
Info << "Creating the Boussinesq buoyancy model" << endl;
vector hRef_ = vector::zero;
label pRefCellGlobal = 0;
label pRefCellProc = 0;
if (pRefCell != -1)
{
    hRef_ = mesh.C()[pRefCell];
    pRefCellGlobal = pRefCell;
    pRefCellProc = Pstream::myProcNo();
}
reduce(hRef_,sumOp<vector>());
reduce(pRefCellGlobal,sumOp<label>());
reduce(pRefCellProc,sumOp<label>());
dimensionedVector hRef("hRef",dimLength,hRef_);
Info << "Pressure reference location is " << hRef.value() << " m at cell " << pRefCellGlobal << " on processor " << pRefCellProc << endl;
    
buoyancyModel Boussinesq(T, TRef, hRef);


// Adjust the background pressure level to agree with the reference pressure specified
// at a certain location.
Info << "Adjusting the background pressure level to agree with reference level" << endl;
#include "adjustPressureLevel.H"


// The quantity p_rgh is required on cell faces (flux) to project the velocity flux given the
// updated pressure field, so add it to the list kept internally in OpenFOAM.
mesh.setFluxRequired(p_rgh.name());


// Allows the user to set fvOptions, such as additional source terms.
#include "createFvOptions.H"


// Allows the user to apply body force with Rhie-Chow type interpolation for stabilization. See:
//
//    Rethore, P-E. M., & Sørensen, N. N., "Actuator Disc Model Using a Modified Rhie-Chow/SIMPLE Pressure
//    Correction Algorithm. Comparison with Analytical Solutions." In Conference proceedings (online) European
//    Wind Energy Association (EWEA), 2008.
//
bool bodyForceStabilization = mesh.schemesDict().lookupOrDefault<bool>("bodyForceStabilization",true);


// Create an object of the horizontalWindTurbineArray class if there is to be a turbine array
turbineModels::horizontalAxisWindTurbinesALMOpenFAST turbines(U);
