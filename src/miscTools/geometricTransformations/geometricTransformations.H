/*---------------------------------------------------------------------------*\
This file was modified or created at the National Renewable Energy
Laboratory (NREL) on January 6, 2012 in creating the SOWFA (Simulator for
Offshore Wind Farm Applications) package of wind plant modeling tools that
are based on the OpenFOAM software. Access to and use of SOWFA imposes
obligations on the user, as set forth in the NWTC Design Codes DATA USE
DISCLAIMER AGREEMENT that can be found at
<http://wind.nrel.gov/designcodes/disclaimer.html>.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Functions to do rotations and coordinate system transformations.

SourceFiles
    geometricTransformations.C

\*---------------------------------------------------------------------------*/

#ifndef geometricTransformations_H
#define geometricTransformations_H

#include "scalar.H"
#include "primitiveFieldsFwd.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    vector rotateVector
    (
        vector v, 
        vector translation, 
        vector axis, 
        scalar angle
    );
    
    vector transformGlobalCartToLocalCart
    (
        vector v, 
        vector xP, 
        vector yP, 
        vector zP
    );

    vectorField transformGlobalCartToLocalCart
    (
        vectorField v,
        vector xP,
        vector yP,
        vector zP
    );
    
    vector transformLocalCartToGlobalCart
    (
        vector vP, 
        vector xP, 
        vector yP, 
        vector zP
    );

    vectorField transformLocalCartToGlobalCart
    (
        vectorField vP,
        vector xP,
        vector yP,
        vector zP
    );
    
    vector transformCartToCyl(vector v);
    
    vector transformCylToCart(vector c);
    
    vector transformGlobalCartToRotorLocalCart
    (
        vector v, 
        vector rotorOrigin,
        vector rotorAxis
    );
    
    vector transformGlobalCartToRotorLocalCyl
    (
        vector v,
        vector rotorOrigin,
        vector rotorAxis
    );
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
