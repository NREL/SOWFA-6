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
    Foam::hurricaneSource

Description
    Source terms use to simulate hurricane flow away from the eyewall following
    G. Bryan's "simple" method.

SourceFiles
    hurricaneSource.C

\*---------------------------------------------------------------------------*/

#ifndef hurricaneSource_H
#define hurricaneSource_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "meshPlanes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class hurricaneSource Declaration
\*---------------------------------------------------------------------------*/

class hurricaneSource
{
    // Private data

        //- Constants
            //- Runtime pointer
            const Time& runTime_;
    
            //- Mesh pointer
            const fvMesh& mesh_;

            //- Velocity field pointer
            const volVectorField& U_;

        //- Dictionary
        const IOdictionary& dict_;

        //- Vertical mesh planes for horizontal averaging
        meshPlanes zPlanes_;
        
        //- Body force field of terms specific to hurricane.
        volVectorField source_;

        bool active_;

        //- Variables related to obtaining the planar-averaged velocity.
        label nLevels;
        List<List<label> > cellsInPlane;
        List<scalar> planeHeights;
        List<vector> Ubar;

        //- Parameters for hurricane source terms.
        scalar R;
        scalar V_surface;
        scalar dVdR_surface;
        scalar dVdz;
        scalar dVdRdz;




    // Private Member Functions

        //- Read a subdict
        void read();

        void setupPlanarAveraging();


public:

    //- Declare name of the class and its debug switch
    ClassName("hurricaneSource");


    //- Constructors
    hurricaneSource
    (
        const IOdictionary& dict,
        const volVectorField& U
    );


    //- Destructor
    virtual ~hurricaneSource();


    //- Update sponge layer total force
    void update();

    //- Return force
    volVectorField& source()
    {
        return source_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

