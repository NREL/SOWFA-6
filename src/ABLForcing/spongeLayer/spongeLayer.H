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
    Foam::spongeLayer

Description
    Sponge layer for damping out vertically propagating gravity waves
    Possible sponge layer types are "Rayleigh", "viscous" or "none".
    - The "Rayleigh" type means that the damping term is computed as nu*(u_ref-u)
      The viscosity coefficient nu has dimensions of 1/s
    - The "viscous" type means that the damping term is computed as nu * Lapl(u)
      The viscosity coefficient nu has dimensions of m**2/s
    - The "none" type means no damping is added

SourceFiles
    spongeLayer.C

\*---------------------------------------------------------------------------*/

#ifndef spongeLayer_H
#define spongeLayer_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "wallDist.H"
#include "interpolateXY.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class spongeLayer Declaration
\*---------------------------------------------------------------------------*/

class spongeLayer
{
    // Private data

        //- Constants
            //- Name
            const word name_;

            //- Runtime pointer
            const Time& runTime_;
    
            //- Mesh pointer
            const fvMesh& mesh_;
    
            //- Velocity field pointer
            const volVectorField& U_;

        //- Dictionary
        dictionary spongeDict_;
        
        //- List of entries in the subdict
        wordList spongesList_;

        //- Reference velocity field
        uniformDimensionedVectorField Uref_;

        //- Viscosity and time scale fields
        volScalarField viscosity_;
        volScalarField viscosityH_;
        volScalarField viscosityV_;
        volScalarField tauH_;
        volScalarField tauV_;
        volScalarField currentViscosity_;

        //- Body force field of the sponge layer
        volVectorField bodyForce_;

        //- Sponge layer parameters
        // Type of damping -- `Rayleigh` or `viscous`
        word type_;
        
        // Patch in which the damping will be applied
        // If specifying patch, specify width as well and the layer start and 
        // end location will be determined automatically. If using standard
        // patch names (west, east, north, south, upper), neither coordIndex
        // nor direction are necessary to be specified. If not specifying the
        // patch (i.e. a non-conventional layer), then all the other parameters
        // need to be speficied
        word patch_;

        // Width of the layer
        scalar width_;

        // Name of the current sponge (for printing purposes)
        word currentSponge_;

        // Parameters of a generic layer. These do not need to be specified if
        // a patch is given. If a non-conventional layer is needed, then:
        // startLocation_: coordinate of the start of the layer
        // coordIndex_: 0, 1, or 2, indicating the direction (x, y, z)
        // direction_: `stepUp` or `stepDown`, indicating direction of cosine
        scalar startLocation_;
        label coordIndex_;
        word direction_;

        // Max damping coefficient
        scalar dampCoeffMax_;
        
        // Components of the velocity to be damped -- `horizontal` or `vertical`
        word dampingComp_;

        // Velocity in which horizontal damping will be damped towards
        scalar Ux_;
        scalar Uy_;

        // Fraction of the width in which the cosine distribution will be
        // applied. The rest of the layer will contain a constant value of 
        // the damping coefficient, given by DampCoeffMax
        // E.g.: if cosFraction = 0.8, then the cosine distribution will
        // take place in 80% of the layer width and the ramaining outer 20% 
        // of layer width will have the maximum damping coefficient
        scalar cosFraction_; 

        // Vertical filter flat
        bool verticalFilter_;

        // Start location of the vertical filter (abs or agl)
        scalar vertFiltStartHeight_;

        // Thickness of the smooth transition of the vertical filter
        scalar vertFiltCosThickness_;

        // Whether to use abs or agl height for vertical filters
        bool useWallDistZ_;

        // Bool to indicate if layer has time-varying properties and thus
        // needs to be updated at every time step
        bool needsUpdating_; 

    // Private Member Functions

        //- Read a subdict
        void readSingleSpongeSubdict_(int s);

        //- Loop through all sponges getting updated viscosity field location
        void getViscosityField_();

        //- Zero out all tau/viscosity fields for updated locations
        void clearViscosityFields_();

        //- Determine viscosity of current layer
        void calculateCurrentSpongeViscosity_();

        //- Add sponge to the viscosity field
        void addSponge_();

        // - Apply vertical filter to the side sponges
        void applyVerticalFilter_();

        // - Fix viscosity field for overlapping regions
        void adjustOverlappingVisc_();

        // - Read simple <time> <value> table and interpolate in time
        scalar readTableAndInterpolate_(dictionary& dict, word parameterTableName);

public:

    //- Declare name of the class and its debug switch
    ClassName("spongeLayer");


    //- Constructors
    spongeLayer
    (
        const word& name,
        const volVectorField& U
    );


    //- Destructor
    virtual ~spongeLayer();


    // Public Member functions

        //- Update sponge layer total force
        void update();

        //- Return force
        volVectorField& force()
        {
            return bodyForce_;
        }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

