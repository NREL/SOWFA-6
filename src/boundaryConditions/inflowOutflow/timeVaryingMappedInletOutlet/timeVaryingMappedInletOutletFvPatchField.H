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
    Foam::timeVaryingMappedInletOutletFvPatchField

Description
    This boundary condition is like inletOutlet in which if the flow is
    outflow, it is set to zero gradient; if the flow is inflow, the value is
    set.  Rather than set it to a single value as with inletOutlet, the value
    is set as timeVaryingMappedFixedValue would set the value using a saved
    time history of spatially varying data.  This is of use in atmospheric
    simulations in which some boundaries may have inflow and outflow faces, and
    the inflow comes from pre-generated, saved turbulent flow.

    As with timeVaryingMappedFixedValue, on inflow faces, this boundary
    condition interpolates the values from a set of supplied points in space
    and time.  Supplied data should be specified in the
    constant/boundaryData/\<patch name\>/ directory where::
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    Alternatively the names and locations of the points and field files may be
    specified explicitly via the optional dictionary entries:
      - dataDir \<optional top-level directory of the points and field data>;
      - points \<optional path including name of points file relative to
                 dataDir\>;
      - sample \<optional name of the sub-directory in the time directories
                containing the fields\>;
    This is particularly useful when mapping data from another case for which
    the \c sample \c functionObject is used to obtain the patch field data for
    mapping.

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first three points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    Inflow values are interpolated linearly between times.

Usage
    \table
        Property     | Description                | Required    | Default value
        setAverage   | Switch to activate setting of average value | no | false
        perturb      | Perturb points for regular geometries | no | 1e-5
        fieldTable   | Alternative field name to sample | no| this field name
        mapMethod    | Type of mapping            | no | planarInterpolation
        offset       | Offset to mapped values    | no | Zero
        dataDir      | Top-level directory of the points and field data \\
                         | no | constant/boundaryData/\<patch name\>
        points       | Path including name of points file relative to dataDir \\
                         | no | points
        sample       | Name of the sub-directory in the time directories \\
                         containing the fields | no | ""
    \endtable

    \verbatim
    <patch name>
    {
        type            timeVaryingMappedInletOutlet;
    }
    \endverbatim

    \verbatim
    <patch name>
    {
        type            timeVaryingMappedInletOutlet;
        dataDir         "../<source case name>/postProcessing/sample";
        points          "0/<sample name>/faceCentres";
        sample          <sample name>;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::inletOutletFvPatchField
    Foam::Function1Types

SourceFiles
    timeVaryingMappedInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedInletOutletFvPatchField_H
#define timeVaryingMappedInletOutletFvPatchField_H

#include "mixedFvPatchFields.H"
#include "surfaceFields.H" // to properly call lookupPatchField
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class timeVaryingMappedInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedInletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;

    // Private data

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- Optional name of the data directory,
        //  defaults to "constant/boundaryData/<patch name>"
        fileName dataDir_;

        //- Optional name of the points file, defaults to "points"
        fileName pointsName_;

        //- Optional name of the sample sub-directory from which the field are
        //  read, defaults to ""
        fileName sampleName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- Does this boundary fix values so that they cannot be adjusted by adjustPhi
        Switch fixesValue_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated values from startSampleTime
        Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        Type startAverage_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated values from endSampleTime
        Field<Type> endSampledValues_;

        //- If setAverage: end average value
        Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;

        //- Helper function to find the field files
        fileName findFieldFile(const word& timeName) const;

        //- Find boundary data in between current time and interpolate
        void checkTable();


public:

    //- Runtime type information
    TypeName("timeVaryingMappedInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedInletOutletFvPatchField
        //  onto a new patch
        timeVaryingMappedInletOutletFvPatchField
        (
            const timeVaryingMappedInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedInletOutletFvPatchField
        (
            const timeVaryingMappedInletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedInletOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedInletOutletFvPatchField
        (
            const timeVaryingMappedInletOutletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return startSampledValues
            const Field<Type> startSampledValues()
            {
                 return startSampledValues_;
            }

            //- Returns true if this patch fixes a value.  This is used for
            //  setting a reference value for a Poisson solve (pressure),
            //  and for determining if adjustPhi affects the flux here.
            virtual bool fixesValue() const
            {
                return fixesValue_;
            }


        // Setters

            void setFixesValue(bool b)
            {
                fixesValue_ = b;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingMappedInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
