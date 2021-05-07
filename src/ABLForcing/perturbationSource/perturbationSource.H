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
    Foam::perturbationSource<Type>

Description
    Add perturbations to the flow field (e.g., to initiate turbulence given
    non-turbulent inflow) with semi-random source terms.

Usage

SourceFiles
    perturbationSource.C

\*---------------------------------------------------------------------------*/

#ifndef perturbationSource_H
#define perturbationSource_H

#include "fvCFD.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class perturbationSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class perturbationSource
{
   // private data:

        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

        //- Constants

            //- Name
            const word name_;

            //- Runtime pointer
            const Time& runTime_;

            //- Mesh pointer
            const fvMesh& mesh_;

            //- Field pointer
            const volFieldType& field_;


        //- Perturbation source field
        volFieldType source_;

public:

    //- Runtime type information
    ClassName("perturbationSource");


    // Constructors
    perturbationSource 
    (
        const IOdictionary& dict,
        const word& name,
        const volFieldType& field
    );


    // Destructor
    virtual ~perturbationSource();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
