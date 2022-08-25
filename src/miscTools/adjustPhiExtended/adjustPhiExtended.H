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

InNamespace
    Foam

Description
    For cases which do no have a pressure boundary adjust the balance of
    fluxes to obey continuity.  Return true if the domain is closed. This
    is an extension of the original OpenFOAM-standard adjustPhiExtendedExtended, 
    which ignores cyclic boundaries.  For the concurrent precursor method 
    that feeds precursor flow to the downstream section via a triple-cyclic 
    boundary set up, one of the precursor cyclics should be considered inflow.

SourceFiles
    adjustPhiExtended.C

\*---------------------------------------------------------------------------*/

#ifndef adjustPhiExtended_H
#define adjustPhiExtended_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Adjust the balance of fluxes to obey continuity.
//  For cases which do not have a pressure boundary.
//  Return true if the domain is closed.
bool adjustPhiExtended
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p,
    List<word> additionalInflowPatch
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
