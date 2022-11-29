// Calculate divergence of velocity flux and display
{

    // Compute divergence cell by cell and report statistics.
    volScalarField divPhi = fvc::div(phi);
    volScalarField divPhiMag  = pow(pow(divPhi,2),0.5);
    scalar minLocalPhiContErr = min(divPhiMag).value();
    reduce(minLocalPhiContErr, minOp<scalar>());
    scalar maxLocalPhiContErr = max(divPhiMag).value();
    reduce(maxLocalPhiContErr, maxOp<scalar>());
    scalar avgLocalPhiContErr = divPhiMag.weightedAverage(mesh.V()).value();
    Info << "Continuity Report:" << endl;
    Info << "   -Local Cell Continuity Error:" << endl;
    Info << "       minimum: " << minLocalPhiContErr << endl;
    Info << "       maximum: " << maxLocalPhiContErr << endl;
    Info << "       weighted mean: " << avgLocalPhiContErr << endl;


    // Compute global divergence over domain boundaries and report statistics.
    Info << "   -Boundary Flux:" << endl;


    List<bool> isCoupledInflow(phi.boundaryField().size(),false);
    forAll(phi.boundaryField(), patchi)
    {
        word patchName = phi.mesh().boundary()[patchi].name();
        forAll(coupledInflowPatches,j)
        {
            if (patchName == coupledInflowPatches[j])
            {
                isCoupledInflow[patchi] = true;
            }
        }
    }

    scalar globalSumPhiBoundary = 0.0;
    scalar globalSumAreaBoundary = 0.0;
    forAll(phi.boundaryField(), patchi)
    {
        // Compute flux on all boundaries except processor boundaries (this is
        // updated from previously where flux was computed for all uncoupled
        // boundaries, i.e., coupled=cyclic and processor, as commented out 
        // below).
      //if ( !mesh.boundary()[patchi].coupled() || isCoupledInflow[patchi] )
        if ( !isA<processorFvPatch>(mesh.boundary()[patchi]) || isCoupledInflow[patchi] )
        {
            word boundaryName = mesh.boundary()[patchi].name();
            scalar sumPhiBoundary = 0.0;
            scalar sumAreaBoundary = 0.0;
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
            forAll(phip,i)
            {
                sumPhiBoundary += phip[i];
                sumAreaBoundary += mesh.boundary()[patchi].magSf()[i];
            }

            reduce(sumPhiBoundary, sumOp<scalar>());
            reduce(sumAreaBoundary, sumOp<scalar>());

            globalSumPhiBoundary += sumPhiBoundary;
            globalSumAreaBoundary += sumAreaBoundary;

            Info << "       " << boundaryName << " - flux: " << sumPhiBoundary << tab << "/ area: " << sumAreaBoundary << endl;
        }
    }

    Info << "       total - flux: " << globalSumPhiBoundary << tab << "/ area: " << globalSumAreaBoundary << nl << endl;
}
