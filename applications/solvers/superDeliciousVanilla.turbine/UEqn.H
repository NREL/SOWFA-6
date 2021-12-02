    // Solve the momentum equation
    
    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    Boussinesq.updateBuoyancyTerm();

    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      + fvc::div(Rwall)
      - momentumSourceTerm.force()
      - Coriolis.force()
      - spongeLayers.force()
      - (1.0 - scalar(bodyForceStabilization)) * turbines.force()
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

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
                  + Boussinesq.buoyancyTerm()
                  + scalar(bodyForceStabilization) * fvc::flux(turbines.force())/mesh.magSf()
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
