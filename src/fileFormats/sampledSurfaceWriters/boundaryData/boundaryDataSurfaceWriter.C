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

#include "boundaryDataSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "Time.H"
#include "pointIOField.H"
#include "primitivePatch.H"

#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(boundaryDataSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::boundaryDataSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    // geometry: rootdir/surfaceName/"points"
    // field:    rootdir/surfaceName/time/field


    const fileName baseDir(outputDir.path()/surfaceName);
    const fileName timeName(outputDir.name());


    mkDir(baseDir);
    OFstream os(baseDir/"points");


    if (isNodeValues)
    {
        if (verbose)
        {
            Info << "Writing points to " << baseDir/"points" << endl;
        }

        os << points;
    }
    else
    {
        if (verbose)
        {
            Info << "Writing face centres to " << baseDir/"points" << endl;
        }

        pointField faceCentres(faces.size(),point::zero);

        forAll(faces, facei)
        {
            faceCentres[facei] = faces[facei].centre(points);
        }

        os << faceCentres;
    }


    // Write field
    {
          fileName valsDir(baseDir/timeName);
          mkDir(valsDir);
          OFstream os(valsDir/fieldName);
          os << values;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryDataSurfaceWriter::boundaryDataSurfaceWriter()
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryDataSurfaceWriter::~boundaryDataSurfaceWriter()
{}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryDataSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    // geometry: rootdir/surfaceName/"points"
    // field:    rootdir/surfaceName/time/field

    const fileName baseDir(outputDir.path()/surfaceName);
    const fileName mainDir(outputDir.path());
    const fileName timeName(outputDir.name());


    rmDir(mainDir/timeName);


    // Write points
    mkDir(baseDir);
    OFstream os(baseDir/"points");

    if (verbose)
    {
        Info << "Writing points to " << baseDir/"points" << endl;
    }

    os << points;
}

// Field writing methods
defineSurfaceWriterWriteFields(Foam::boundaryDataSurfaceWriter);

// ************************************************************************* //
