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
    
    if (subDict.empty())
    {
        Info << "No subdict..." << endl;
        
    }

    if (!subDict.found("V"))
    {
        FatalErrorInFunction << "Must specify hurricane characteristic velocity, 'V', at this location."
                                 << abort(FatalError);
    }
    
    if (!subDict.found("R"))
    {
        FatalErrorInFunction << "Must specify radius from center of hurricane, 'R'."
                                 << abort(FatalError);
    }

    if (!subDict.found("dVdR"))
    {
        FatalErrorInFunction << "Must specify radial characteristic velocity gradient, 'dVdR', at this location"
                                 << abort(FatalError);
    }

    V = subDict.lookupOrDefault<scalar>("V",10.0);
    R = subDict.lookupOrDefault<scalar>("R",40.0E3);
    dVdR = subDict.lookupOrDefault<scalar>("dVdR",0.0);
}


void Foam::hurricaneSource::setupPlanarAveraging()
{
    // Initialize the planar averaging.
    zPlanes_.initialize();

    nLevels = zPlanes_.numberOfPlanes();

    cellsInPlane = zPlanes_.planesCellList();
}


void Foam::hurricaneSource::update()
{
    // Get planar average velocity at each height.
    Ubar = zPlanes_.average(U_);

    // Get the Coriolis f-factor.
    const scalar f = Coriolis_.f();

    // Update the source term.
    for (int i = 0; i < nLevels; i++)
    {
        forAll(cellsInPlane[i], j)
        {
            source_[j].x() =  Foam::sqr(Ubar[i].x()) / R  
                             +Ubar[i].y() * (V/R)
                             -(f*V + (Foam::sqr(V) / R));
            source_[j].y() = -Ubar[i].x() * dVdR 
                             -Ubar[i].x() * (V/R);
            source_[j].z() = 0.0;
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hurricaneSource::hurricaneSource
(
    const IOdictionary& dict,
    const volVectorField& U,
    const CoriolisForce& Coriolis
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    Coriolis_(Coriolis),

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
    )
{

read();
setupPlanarAveraging();
update();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hurricaneSource::~hurricaneSource()
{}


// ************************************************************************* //
