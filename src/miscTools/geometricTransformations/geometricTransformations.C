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

\*---------------------------------------------------------------------------*/

#include "geometricTransformations.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::rotateVector
(
    vector v, 
    vector origin, 
    vector axis,
    scalar angle
)
{
    // Rotate a point about an origin and about an arbitrary axis by a given angle.

    // Declare and define the rotation matrix.
    tensor RM;
    RM.xx() = Foam::sqr(axis.x()) + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle); 
    RM.xy() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle); 
    RM.xz() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y() * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle); 
    RM.yy() = Foam::sqr(axis.y()) + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z() * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z() * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z()) + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract rotation point
    // off the point to be rotated.
    v = v - origin;

    // Perform the rotation.
    v = RM & v;

    // Return the rotated point to its new location relative to the rotation point.
    v = v + origin;

    return v;
}



Foam::vector Foam::transformGlobalCartToLocalCart
(
    vector v, 
    vector xP, 
    vector yP, 
    vector zP
)
{
    // Transform from the global Cartesian (x,y,z) system into the local 
    // Cartesian (x',y',z') system using v' = T'v
    //
    //    x' is aligned with the flow
    //    y' is the cross product of z' and x'
    //    z' is in the boundary face normal direction
    //
    // These vectors are unit vectors.  The vectors make up the rows of
    // the rotation matrix T', which rotates from (x,y,z) to (x',y',z').

    scalar zPMag;
    zPMag = mag(zP);
    zP = zP/zPMag;

    scalar xPMag;
    xPMag = mag(xP);
    xP = xP/xPMag;

    scalar yPMag;
    yPMag = mag(yP);
    yP = yP/yPMag;

    // Create T'
    tensor TP;
    TP.xx() = xP.x();
    TP.xy() = xP.y();
    TP.xz() = xP.z();
    TP.yx() = yP.x();
    TP.yy() = yP.y();
    TP.yz() = yP.z();
    TP.zx() = zP.x();
    TP.zy() = zP.y();
    TP.zz() = zP.z();

    // Transform the vector from Cartesian to local
    vector vP = TP & v;

    return vP;
}



Foam::vectorField Foam::transformGlobalCartToLocalCart
(
    vectorField v,
    vector xP,
    vector yP,
    vector zP
)
{
    vectorField vP(v);
    forAll(vP,i)
    {
        vP[i] = transformGlobalCartToLocalCart(v[i],xP,yP,zP);
    }

    return vP;
}



Foam::vector Foam::transformLocalCartToGlobalCart
(
    vector vP,
    vector xP, 
    vector yP, 
    vector zP
)
{
    // Transform from the local Cartesian (x',y',z') system to the global 
    // Cartesian (x,y,z) system using v = Tv'
    //
    //    x' is aligned with the flow
    //    y' is the cross product of z' and x'
    //    z' is in the boundary face normal direction
    //
    // These vectors are unit vectors.  The vectors make up the rows of
    // the rotation matrix T', which rotates from (x,y,z) to (x',y',z').
    // T can be recovered from T' because it is the inverse. T' is
    // such that the (T')^-1 = transpose(T') because it is made up of
    // orthogonal basis vectors. 

    scalar zPMag;
    zPMag = mag(zP);
    zP = zP/zPMag;

    scalar xPMag;
    xPMag = mag(xP);
    xP = xP/xPMag;

    scalar yPMag;
    yPMag = mag(yP);
    yP = yP/yPMag;

    // Create T'
    tensor TP;
    TP.xx() = xP.x();
    TP.xy() = xP.y();
    TP.xz() = xP.z();
    TP.yx() = yP.x();
    TP.yy() = yP.y();
    TP.yz() = yP.z();
    TP.zx() = zP.x();
    TP.zy() = zP.y();
    TP.zz() = zP.z();

    // Create T
    tensor T = TP.T();

    // Transform the vector from Cartesian to local
    vector v = T & vP;

    return v;
}



Foam::vectorField Foam::transformLocalCartToGlobalCart
(
    vectorField vP,
    vector xP,
    vector yP,
    vector zP
)
{
    vectorField v(vP);
    forAll(v,i)
    {
        v[i] = transformLocalCartToGlobalCart(vP[i],xP,yP,zP);
    }
   
    return v;
}



Foam::vector Foam::transformCartToCyl(vector v)
{
    // Transform from the Cartesian (x,y,z) system to the cylindrical 
    // (r,theta,x) system.

    // Convert to cylindrical coordinates.
    vector c = vector::zero;
    c.x() = Foam::sqrt(Foam::sqr(v.y()) + Foam::sqr(v.z()));
    c.y() = Foam::atan2(v.z(),v.y());
    c.z() = v.x();

    return c;
}



Foam::vector Foam::transformCylToCart(vector c)
{
    // Transform from the cylindrical (r,theta,x) system to the Cartesian
    // (x,y,z) system.

    // Convert to Cartesian coordinates.
    vector v = vector::zero;
    v.x() = c.z();
    v.y() = c.x()*Foam::cos(c.y());
    v.z() = c.x()*Foam::sin(c.y());

    return v;
}



Foam::vector Foam::transformGlobalCartToRotorLocalCart
(
    vector v, 
    vector rotorOrigin, 
    vector rotorAxis
)
{
    // Transform a point from global Cartesian coordinates to rotor local Cartesian coordinates.

    // Get the vector relative to the rotor apex.
    vector pCartGlobal = v - rotorOrigin;

    // Define the orientation of the local Cartesian coordinate system such that x is along the shaft,
    // z is up, but normal to the shaft, and y is orthogonal to x and z.
    vector xP = rotorAxis;
    xP /= mag(xP);

    vector zP = vector::zero;
    zP.z() = 1.0;

    vector yP = -(xP ^ zP);
    yP /= mag(yP);

    zP = xP ^ yP;

    // Transform to the local Cartesian system.
    vector pCartLocal = transformGlobalCartToLocalCart(pCartGlobal,xP,yP,zP);

    return pCartLocal;
}



Foam::vector Foam::transformGlobalCartToRotorLocalCyl
(
    vector v,
    vector rotorOrigin,
    vector rotorAxis
)
{
    // Transform a point from global Cartesian coordinates to rotor local cyclindrical coordinates.
    
    // Transform from global Cartesian to rotor local Cartesian.
    vector pCartLocal = transformGlobalCartToRotorLocalCart(v,rotorOrigin,rotorAxis);

    // Transform to the local cylindrical system.
    vector pCylLocal = transformCartToCyl(pCartLocal);

    return pCylLocal;
}
    
// ************************************************************************* //

