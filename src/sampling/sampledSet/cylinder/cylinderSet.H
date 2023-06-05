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

Class
    Foam::cylinderSet

Description
    Samples a regular array of points on a cylinder.

Usage
    \table
        Property           | Description                             | Req'd? | Default
        axisPointStart     | The starting point of the cylinder axis | yes    |
        axisPointEnd       | The ending point of the cylinder axis   | yes    |
        radiusStart        | The starting radius of the cylinder     | yes    |
        radiusEnd          | The ending radius of the cylinder       | yes    |
        thetaStart         | The starting azimuth (in degrees)       | yes    |
        thetaEnd           | The ending azimuth (in degrees)         | yes    |
        pointsDensity      | The number of points in the box         | yes    |
        axis               | The coordinate axis that is written     | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type               cylinder;
        axisPointStart    (-15.0 0.0 0.0);
        axisPointEnd      ( 15.0 0.0 0.0);
        radiusStart        25.0;
        radiusEnd          50.0;
        thetaStart         0.0;
        thetaEnd          270.0;
        pointsDensity     (31 7 26);
        axis               xyz;
    }
    \endverbatim

SourceFiles
    cylinderSet.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderSet_H
#define cylinderSet_H

#include "sampledSet.H"
#include "labelVector.H"
#include "DynamicList.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                           Class cylinderSet Declaration
\*---------------------------------------------------------------------------*/

class cylinderSet
:
    public sampledSet
{
    // Private data

        //- Origin (x, y, z) in global cartesian co-ordinates
        point axisPointStart_;
        point axisPointEnd_;

        //- Radius of cylinder
        scalar radiusStart_;
        scalar radiusEnd_;

        //- Azimuth start and end locations in degrees and radians.
        scalar thetaStartDeg_;
        scalar thetaEndDeg_;

        //- Number of points in each direction.
        labelVector pointsDensity_;



    // Private Member Functions

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("cylinder");


    // Constructors

        //- Construct from components
        cylinderSet
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
        );

        //- Construct from dictionary
        cylinderSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cylinderSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
