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
    Foam::boundaryDataSurfaceWriter

Description
    A surfaceWriter for outputting to a form useable for the
    timeVaryingMapped boundary condition. This reads the data from
    constant/boundaryData/\<patch\> directory

    Typical way of working:
    - use a sampledSurface of type 'patch' (to sample a patch):
    \verbatim
    surfaces
    {
        type            surfaces;
        surfaceFormat   boundaryData;
        fields          ( p );
        surfaces
        (
            outlet
            {
                type            patch;
                patches         (outlet);
                interpolate     false;
            }
        );
    }
    \endverbatim

    - write using this writer.
    - move postProcessing/surfaces/outlet to constant/boundaryData/outlet
      in your destination case.
    - use a timeVaryingMappedFixedValue condition to read and interpolate
      the profile:
        type            timeVaryingMappedFixedValue;
        setAverage      false;  // do not use read average
        offset          0;      // do not apply offset to values

    Note:
    - with 'interpolate false' the data is on the face centres of the
      patch. Take care that a 2D geometry will only have a single row
      of face centres so might not provide a valid triangulation
      (this is what timeVaryingMappedFixedValue uses to do interpolation)
      (Alternatively use timeVaryingMappedFixedValue with mapMethod 'nearest')

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    where the geometry is written as:
    \verbatim
    rootdir
    `-- surfaceName
        `-- "points"
    \endverbatim

    and field data:
    \verbatim
    rootdir
    `-- surfaceName
        |-- "points"
        `-- timeName
            `-- field
    \endverbatim

SourceFiles
    boundaryDataSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryDataSurfaceWriter_H
#define boundaryDataSurfaceWriter_H

#include "surfaceWriter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class boundaryDataSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class boundaryDataSurfaceWriter
:
    public surfaceWriter
{
    // Private Member Functions

        //- Templated write operation
        template<class Type>
        void writeTemplate
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;

public:

    //- Runtime type information
    TypeName("boundaryData");


    // Constructors

        //- Construct null
        boundaryDataSurfaceWriter();


    //- Destructor
    virtual ~boundaryDataSurfaceWriter();


    // Member Functions

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry()
        {
            return true;
        }

        //- Write single surface geometry to file.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const;


        //- Write scalarField for a single surface to file.
        //  One value per face or vertex.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<scalar>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<vector>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<sphericalTensor>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<symmTensor>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex.
        virtual void write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const pointField& points,
            const faceList& faces,
            const word& fieldName,
            const Field<tensor>& values,
            const bool isNodeValues = false,
            const bool verbose = false
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
