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

#include "cylinderSet.H"
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
    defineTypeNameAndDebug(cylinderSet, 0);
    addToRunTimeSelectionTable(sampledSet, cylinderSet, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::cylinderSet::calcSamples
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
  
    // The cylinder may be arbitrarily oriented in x-y space. Get necessary orientation
    // information.  (Need to update for true fully arbitrary orientation.)  
    //  - this is the axis down the length of the cylinder.
    vector cylinderAxis = axisPointEnd_ - axisPointStart_;

    //  - this is the length of the cylinder.
    scalar cylinderLength = mag(cylinderAxis);

    //  - make the cyliner axis a unit vector.
    cylinderAxis /= mag(cylinderAxis);


    //  - this is an axis of a radial of the cylinder.  Pointed from the cylinder axis
    //  - outward along the radius.
    vector radialAxis = Zero;
    radialAxis.x() = -cylinderAxis.y();
    radialAxis.y() =  cylinderAxis.x();
    radialAxis /= mag(radialAxis);

    // Convert theta to radians.
    scalar thetaStartRad_ = thetaStartDeg_ * (Foam::constant::mathematical::pi/180.0);
    scalar thetaEndRad_ = thetaEndDeg_ * (Foam::constant::mathematical::pi/180.0);

    // Compute the sampling point spacing.
    const scalar dx = (pointsDensity_.x() > 1) ? cylinderLength/(pointsDensity_.x() - 1) : 0;
    const scalar dt = (pointsDensity_.y() > 1) ? (thetaEndRad_ - thetaStartRad_)/(pointsDensity_.y() - 1) : 0;
    const scalar dtDeg = dt*(180.0/Foam::constant::mathematical::pi);
    const scalar dr = (pointsDensity_.z() > 1) ? (radiusEnd_ - radiusStart_)/(pointsDensity_.z() - 1) : 0;
    scalar avgRes = (dx + dr + dt)/3.0;

    Info << nl << endl;
    Info << "Cylinder sampling:" << endl;
    Info << "  -resolution = (" << dx << " m, " << dtDeg << " deg, " << dr << " m)" << endl;
    Info << "  -cylinderLength = " << cylinderLength << " m" << endl;
    Info << "  -radius range = " << radiusStart_ << " to " << radiusEnd_ << " m" << endl;
    Info << "  -azimuthal range = " << thetaStartDeg_ << " to " << thetaEndDeg_ << " deg" << endl;
    Info << nl << endl;


    // Compute and locate the points.
    label pointCount = 0;

    for (label k = 0; k < pointsDensity_.z(); k++)
    {
        // Get the radius at this k-level of the sampling mesh.
        scalar radius = radiusStart_ + (k*dr);

        for (label j = 0; j < pointsDensity_.y(); j++)
        {

            // Get the azimuth at this j-level of the sampling mesh.
            scalar theta = thetaStartRad_ + (j*dt);

            // Get this particular radial axis.  Start with the basic radial
            // axis and rotate.
            vector r = radialAxis;
            r = r * cos(theta) + (cylinderAxis^r)*sin(theta);
            r /= mag(r);
 
            for (label i = 0; i < pointsDensity_.x(); i++)
            {
                vector origin = axisPointStart_ + (i*dx)*cylinderAxis;
                point pt(origin + radius*r);

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


void Foam::sampledSets::cylinderSet::genSamples()
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

Foam::sampledSets::cylinderSet::cylinderSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& axisPointStart,
    const point& axisPointEnd,
    const scalar& radiusStart,
    const scalar& radiusEnd,
    const scalar& thetaStart,
    const scalar& thetaEnd,
    const labelVector& pointsDensity
)
:
    sampledSet(name, mesh, searchEngine, axis),
    axisPointStart_(axisPointStart),
    axisPointEnd_(axisPointEnd),
    radiusStart_(radiusStart),
    radiusEnd_(radiusEnd),
    thetaStartDeg_(thetaStart),
    thetaEndDeg_(thetaEnd),
    pointsDensity_(pointsDensity)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::sampledSets::cylinderSet::cylinderSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    axisPointStart_(dict.lookup("axisPointStart")),
    axisPointEnd_(dict.lookup("axisPointEnd")),
    radiusStart_(readScalar(dict.lookup("radiusStart"))),
    radiusEnd_(readScalar(dict.lookup("radiusEnd"))),
    thetaStartDeg_(readScalar(dict.lookup("thetaStart"))),
    thetaEndDeg_(readScalar(dict.lookup("thetaEnd"))),
    pointsDensity_(dict.lookup("pointsDensity"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::cylinderSet::~cylinderSet()
{}


// ************************************************************************* //
