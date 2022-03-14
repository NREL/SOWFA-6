    // Solve the momentum equation
    
    // Set up the solution matrix
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    ==
      - fvc::div(Rwall)
      + Coriolis.force()
      + momentumGeoMesoTerm.force()
      + momentumSpongeLayers.force()
      + momentumPerturbationZones.source()
      + fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    {
    volScalarField aa = UEqn.A();
    volVectorField hh = UEqn.H();
    Info << "UEqn.A():" << endl << aa[0] << endl;
    Info << "UEqn.H():" << endl << hh[0] << endl;
    }

/*
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p_rgh)
                  + Boussinesq.source()
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

    if (pimple.finalPimpleIter() && pimple.finalPISOIter())
    {
        momentumPerturbationZones.perturb(U);
    }
*/
