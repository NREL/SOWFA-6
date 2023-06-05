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

#include "arrayStructured.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "transform.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(arrayStructured, 0);
    addToRunTimeSelectionTable(sampledSet, arrayStructured, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::arrayStructured::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Set a random number generator.
    Random randGen(label(0));

    // Compute the sampling points spacing.
    const scalar dx = (pointsDensity_.x() > 1) ? spanBox_.x()/(pointsDensity_.x() - 1) : 0;
    const scalar dy = (pointsDensity_.y() > 1) ? spanBox_.y()/(pointsDensity_.y() - 1) : 0;
    const scalar dz = (pointsDensity_.z() > 1) ? spanBox_.z()/(pointsDensity_.z() - 1) : 0;
    scalar avgRes = (dx + dy + dz)/3.0;

    Info << nl << endl;
    Info << "Cartesian Structured Array Sampling:" << endl;
    Info << "  -resolution = (" << dx << " m, " << dy << " m, " << dz << " m)" << endl;
    Info << "  -box size = " << spanBox_ << " m" << endl;
    Info << nl << endl;

    // Compute and locate the points.
    label pointCount = 0;

    for (label k = 0; k < pointsDensity_.z(); k++)
    {
        for (label j = 0; j < pointsDensity_.y(); j++)
        {
            for (label i = 0; i < pointsDensity_.x(); i++)
            {
		// Local Cartesian
		point pt(i*dx, j*dy, k*dz);

		// Global Cartesian
		pt = csys_.globalPosition(pt);

                // Compute a perturbation to the true point location only for the purposes
                // of searching to eliminate ties for points exactly on cell faces on
                // processor boundaries.
                
                vector ptPerturb = 0.001 * avgRes * (2.0*randGen.sample01<vector>() - vector::one);
                const label celli = searchEngine().findCell(pt+ptPerturb);
              //Info << "ptP(" << i << ", " << j <<", " << k << ") = " << pt+ptPerturb << endl; 
              //Info << "pt(" << i << ", " << j <<", " << k << ") = " << pt << endl; 

                if (celli != -1)
                {
                    samplingPts.append(pt);
                    samplingCells.append(celli);
                    samplingFaces.append(-1);
                    samplingSegments.append(0);
                    samplingCurveDist.append(scalar(pointCount));
                }
                pointCount++;
            }
        }
    }
}


void Foam::sampledSets::arrayStructured::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::arrayStructured::arrayStructured
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const cartesianCS& csys,
    const Vector<label>& pointsDensity,
    const Vector<scalar>& spanBox
)
:
    sampledSet(name, mesh, searchEngine, axis),
    csys_(csys),
    pointsDensity_(pointsDensity),
    spanBox_(spanBox)
{
    genSamples();
}



Foam::sampledSets::arrayStructured::arrayStructured
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    csys_(dict),
    pointsDensity_(dict.lookup("pointsDensity")),
    spanBox_(dict.lookup("spanBox"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::arrayStructured::~arrayStructured()
{}

// ************************************************************************* //
