/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adjustPhiExtended.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::adjustPhiExtended
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p,
    List<word> additionalInflowPatch
)
{
    if (p.needReference())
    {
      //Info << "Calling adjustPhiExtended..." << endl;
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        surfaceScalarField::Boundary& bphi =
            phi.boundaryFieldRef();


        // Check to see if a patch is considered inflow.  The conditions for
        // this to be true are that patch specifies mass flow rate (i.e., 
        // 'fixesValue()' returns true) but is not an inletOutlet,  or the 
        // patch has been specifically labeled as inflow in the 
        // 'additionalInflowPatch' list.  This is useful for cases where
        // there is a concurrent precursor simulation that is laterally
        // cyclic, but we want to count one or more of these cyclic patches
        // as an inflow.
        List<bool> isInflow(bphi.size(),false);
        forAll(bphi, patchi)
        {
            word patchName = phi.mesh().boundary()[patchi].name();
            forAll(additionalInflowPatch,j)
            {
                if (patchName == additionalInflowPatch[j])
                {
                    isInflow[patchi] = true;
                }
            }
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            if (Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up))
            {
                isInflow[patchi] = true;
            }
        }


        // Go through and sum up the mass flow rate into the domain, the
        // fixed mass flow rate out of the domain, and the adjustable
        // mass flow rate out of the domain.  Adjustable mass flow rate
        // happens on inletOutlet faces, for example, but fixed mass flow
        // rate out of the domain happens on fixedValue boundaries.
        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled() || isInflow[patchi])
            {
                if (isInflow[patchi])
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = vSmall + sum(mag(phi)).value();

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        // Compute the correction to the adjustable mass flow rate out such that it
        // balances with the net of fixed in and outflow.
        if
        (
            magAdjustableMassOut > vSmall
         && magAdjustableMassOut/totalFlux > small
        )
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            FatalErrorInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }


        // Apply the correction to the adjustable mass flow rate out faces.
        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled() || isInflow[patchi])
            {
                if (!isInflow[patchi])
                {
                    word patchName = phi.mesh().boundary()[patchi].name();
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        return mag(massIn)/totalFlux < small
            && mag(fixedMassOut)/totalFlux < small
            && mag(adjustableMassOut)/totalFlux < small;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
