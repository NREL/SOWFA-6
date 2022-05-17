{
    // The inverse of the momentum equation "A" operator (the matrix diagonal) and
    // its projection from cell centers to cell faces.  Because it is a diagonal matrix,
    // The inverse is just the diagonal matrix of the reciprocals, hence the "r" in the
    // name "rAU".
    volScalarField rAU("rAU", 1.0/UEqn.A());
  //volScalarField rAU("rAU", 1.0/A);
  //volScalarField rAU("rAU",mesh, dimensionedScalar("rAU",dimensionSet(0,0,1,0,0,0,0),1.0/1.5),extrapolatedCalculatedFvPatchScalarField::typeName);
  //volScalarField rAU("rAU",1.0/A);
  //rAU = 1.0/1.5;
  //rAU = 1.0/UEqn.A();
  //rAU.correctBoundaryConditions();
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    // The momentum equation "H" operator is the off-diagonal times the last known U.
    // In the equations it is always multiplied with inv(A).  This is HbyA.
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
  //volVectorField HbyA(rAU*UEqn.H());

    // This is the Boussinesq buoyancy term multipled by the inverse of the A operator.
    Boussinesq.updateBuoyancyTerm();
    surfaceScalarField phig(rAUf * Boussinesq.source() * mesh.magSf());

    // This is the actuator line/disk body force term multiplied by the inverse of the A operator.
    surfaceScalarField phiF(rAUf * scalar(bodyForceStabilization) * fvc::flux(actuatorTurbineArray.source()));

    // Project HbyA to cell faces and apply a correction for time stepping scheme.
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
    );

    // This part balances global mass flux.  It does it in a temporary field, and then
    // applies the correction indirectly by setting the pressure gradient to be used in
    // the fixedFluxPressure boundary condition on p_rgh, or directly if the zeroGradient
    // boundary condition on p_rgh is used.
    surfaceScalarField phiFixedFlux = 1.0*phiHbyA;
  //surfaceScalarField phigFixedFlux = 1.0*phig;

  //Info << "phig.boundaryFieldRef() before masking: " << endl;
  //forAll(phigFixedFlux.boundaryField(), patchi)
  //{
  //    Info << mesh.boundary()[patchi].name() << endl;
  //    Info << phigFixedFlux.boundaryField()[patchi] << endl;
  //    Info << rAUf.boundaryField()[patchi] << endl;
  //    Info << T.boundaryField()[patchi] << endl;
  //    Info << Boussinesq.rhok().boundaryField()[patchi] << endl;
  //    Info << Boussinesq.source().boundaryField()[patchi] << endl;
  //}
    /*
    label patchFrom = 0;
    label patchTo = 0;
    forAll(phiFixedFlux.boundaryField(), patchi)
    {
        if (mesh.boundary()[patchi].name() == "east.1")
        {
            patchFrom = patchi;
        }
        if (mesh.boundary()[patchi].name() == "east.2")
        {
            patchTo = patchi;
        }
    }
    Info << "patchFrom = " << patchFrom << tab << mesh.boundary()[patchFrom].name() << endl;
    Info << "patchTo   = " << patchTo   << tab << mesh.boundary()[patchTo  ].name() << endl;
    
    forAll(phiFixedFlux.boundaryField()[patchTo],i)
    {
        phiFixedFlux.boundaryFieldRef()[patchTo][i] = phiFixedFlux.boundaryField()[patchFrom][i];
     // else
     // {
     //     phigFixedFlux.boundaryFieldRef()[patchi] = 0.0;
     // }
    }
  //Info << "phig.boundaryFieldRef() after masking: " << endl;
  //forAll(phigFixedFlux.boundaryField(), patchi)
  //{
  //    Info << mesh.boundary()[patchi].name() << endl;
  //    Info << phigFixedFlux.boundaryFieldRef()[patchi] << endl;
  //}
*/

    adjustPhi(phiFixedFlux, U, p_rgh);

    forAll(p_rgh.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p_rgh.boundaryField()[patchi]))
        {
            phiHbyA.boundaryFieldRef()[patchi] = phiFixedFlux.boundaryField()[patchi];
        }
    }

    phiHbyA += phig;
    phiHbyA += phiF;

    /*
    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - phiFixedFlux.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );
    */

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - phiFixedFlux.boundaryField()
      //  + phigFixedFlux.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );


    
    // Non-orthogonal corrector loop.
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        if (activatePressureRefCell)
        {
            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
        }

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

          //pF = rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
          //pF = fvc::reconstruct((phig - p_rghEqn.flux()));
          //pF = HbyA;
          //pF = rAU*UEqn.H();
          //pF = UEqn.H();
          //pF = rAU;

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig + phiF - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "adjustPressureLevel.H"
}
