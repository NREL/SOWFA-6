/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "hurricaneSource.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hurricaneSource, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hurricaneSource::read()
{
    dictionary subDict = dict_.subOrEmptyDict("hurricaneSource");
    
    active_ = subDict.empty() ? false : true;

    if (active_)
    {
        if (!subDict.found("R"))
        {
            FatalErrorInFunction << "Must specify radius from center of hurricane, 'R'."
                                     << abort(FatalError);
        }

        if (!subDict.found("V_surface"))
        {
            FatalErrorInFunction << "Must specify characteristic velocity at the surface, 'V_surface'."
                                     << abort(FatalError);
        }
        
        if (!subDict.found("dVdR_surface"))
        {
            FatalErrorInFunction << "Must specify characteristic velocity gradient with respect to radius at the surface, 'dVdR_surface'."
                                     << abort(FatalError);
        }

        if (!subDict.found("dVdz"))
        {
            FatalErrorInFunction << "Must specify characteristic velocity gradient with height, 'dVdz'."
                                     << abort(FatalError);
        }
        
        if (!subDict.found("dVdRdz"))
        {
            FatalErrorInFunction << "Must specify gradient of characteristic velocity gradient with respect to radius with respect to height, 'dVdRdz'."
                                     << abort(FatalError);
        }

        R = subDict.lookupOrDefault<scalar>("R",40.0E3);
        V_surface = subDict.lookupOrDefault<scalar>("V_surface",40.0);
        dVdR_surface = subDict.lookupOrDefault<scalar>("dVdR_surface",-8.0E-3);
        dVdz = subDict.lookupOrDefault<scalar>("dVdz",-2.0E-3);
        dVdRdz = subDict.lookupOrDefault<scalar>("dVdRdz",1.0E-7);

        Info << "  -using values:" << endl;
        Info << "      R = " << R/1000.0 << " km" << endl;
        Info << "      V at surface = " << V_surface << " m/s" << endl;
        Info << "      dV/dR at surface = " << dVdR_surface << " 1/s" << endl;
        Info << "      dV/dz = " << dVdz << " 1/s" << endl;
        Info << "      (dV/dR)/dz = " << dVdRdz << " 1/(m*s)" << endl;
    }

    else
    {
        Info << "  -no hurricane source terms specified. Skipping..." << endl;   
    }

}


void Foam::hurricaneSource::setupPlanarAveraging()
{
    if (active_)
    {
        // Initialize the planar averaging.
        zPlanes_.initialize();
    
        nLevels = zPlanes_.numberOfPlanes();
    
        cellsInPlane = zPlanes_.planesCellList();
    
        planeHeights = zPlanes_.planeLocationValues();
    }
}


void Foam::hurricaneSource::update()
{
    if (active_)
    {
        // Get planar average velocity at each height.
        Ubar = zPlanes_.average(U_);

        // Update the source term.
        for (int i = 0; i < nLevels; i++)
        {
            vector sourceAtLevel = Zero;

            scalar V = V_surface + (planeHeights[i] * dVdz);
            scalar dVdR = dVdR_surface + (planeHeights[i] * dVdRdz);

            sourceAtLevel.x() =  Foam::sqr(Ubar[i].x()) / R
                                +Ubar[i].y() * (V/R)
                                -(Foam::sqr(V) / R);
            sourceAtLevel.y() = -Ubar[i].x() * dVdR
                                -Ubar[i].x() * (V/R);

          //Info << planeHeights[i] << tab << V_surface << tab << dVdR_surface << tab << V << tab << dVdR << tab << sourceAtLevel << endl;

            forAll(cellsInPlane[i], j)
            {
                source_[cellsInPlane[i][j]] = sourceAtLevel;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hurricaneSource::hurricaneSource
(
    const IOdictionary& dict,
    const volVectorField& U
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    dict_(dict),

    zPlanes_(mesh_),

    // Set the pointer to the velocity field
    U_(U),

    // Initialize the body force field
    source_
    (
        IOobject
        (
            "hurricaneSource",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimVelocity/dimTime,vector::zero)
    ),

    active_(false)
{
    read();
    setupPlanarAveraging();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hurricaneSource::~hurricaneSource()
{}


// ************************************************************************* //
