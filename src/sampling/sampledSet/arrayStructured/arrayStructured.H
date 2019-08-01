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

Class
    Foam::sampledSets::arrayStructured

Description
    Specifieds an ordered array of uniformly distributed sampling points

Usage
    \table
        Property    | Description                            | Req'd? | Default
        box         | The box which contains the samples     | yes    |
        nPoints     | The number of points in each direction | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        arrayStructured;
        box         (0.95 0 0.25) (1.2 0.25 0.5);
        nPoints     (2 4 6);
        axis        x;
    }
    \endverbatim

SourceFiles
    arrayStructured.C

\*---------------------------------------------------------------------------*/

#ifndef arrayStructured_H
#define arrayStructured_H

#include "sampledSet.H"
#include "labelVector.H"
#include "DynamicList.H"
#include "cartesianCS.H"
#include "cylindricalCS.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                           Class arrayStructured Declaration
\*---------------------------------------------------------------------------*/

class arrayStructured
:
    public sampledSet
{
    // Private data

        //- Local Cartesian coordinate system
	cartesianCS csys_;

        //- The number of points across each direction of the array
        labelVector pointsDensity_;

	//- Span box
	vector spanBox_;


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
    TypeName("arrayStructured");


    // Constructors

        //- Construct from components
        arrayStructured
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
	    const word& axis,
	    const cartesianCS& csys,
	    const Vector<label>& pointDensity,
	    const Vector<scalar>& spanBox
        );

	//- Construct from dictionary
	arrayStructured
	(
            const word& name,
	    const polyMesh& mesh,
	    const meshSearch& searchEngine,
	    const dictionary& dict
	);


    //- Destructor
    virtual ~arrayStructured();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
