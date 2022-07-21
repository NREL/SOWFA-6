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
    Foam::perturbationZone<Type>

Description
    Add perturbations to the flow field (e.g., to initiate turbulence given
    non-turbulent inflow) with semi-random source terms.

Usage

SourceFiles
    perturbationZone.C

\*---------------------------------------------------------------------------*/

#ifndef perturbationZone_H
#define perturbationZone_H

#include "scalar.H"
#include "fvCFD.H"
#include "IOdictionary.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class perturbationZone Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class perturbationZone
{
   // private data:

        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

        //- Constants
            const IOdictionary& dict_;

            //- Runtime pointer
            const Time& runTime_;

            //- Mesh pointer
            const fvMesh& mesh_;

            //- Field pointer (field cannot be modified)
            const volFieldType& field_;

            //- Pointer to velocity field
            const volVectorField& U_; 

        //- Field for storing height above the ground
        volScalarField zAgl_;

        //- Random number generator
        Random randGen_;

        //- Perturbation source field
        volFieldType source_;

        //- A perturbation field if the field is directly modified rather
        //  that indirectly modified through source terms.
        volFieldType fieldPerturbation_;

        //- Perturbation subdictionary within the overall dictionary.
        dictionary subDict_;

        //- List of the names of the perturbation zones.
        wordList subDictList_;

        //- Number of perturbation zones.
        label nZones_;

        //- Defines how the zones location is defined, choices are:
        //     -boundary: a zone adjacent to and touching a given boundary.
        //     -arbitraryBox: a box zone that is stand-alone in the domain.
        List<word> locationType_;

        //- For zones where the locationType is "boundary", then use these location
        //  definitions:

            //-- name of boundary to which the zone will be adjacent.
            List<word> adjacentBoundary_;

            //-- thickness of the zone.
            List<scalar> layerThickness_;

            //-- height of the zone.
            List<scalar> layerHeight_;

            //-- use height relative to the ground instead of absolute height
            List<bool> useWallDist_;
            bool useWallDistAny_;

        //- For zones where the locationType is "arbitraryBox", then use these location
        //  definitions:

            //-- origin of the box
            List<vector> boxOrigin_;

            //-- vectors defining the box 
            List<vector> boxVec_i_;
            List<vector> boxVec_j_;
            List<vector> boxVec_k_;

        //- The resolution of the perturbation boxes
            List<vector> res_;

        //- The dimensions of the perturbation boxes
            List<List<label> > dims_;

        //- The mode in which fluctuation magnitude is computed. Option are:
        //  --manual: User specifies max fluctuation magnitude.
        //  --Eckert: Eckert number formulation. 
            List<word> fluctuationMagMode_;

        //- The magnitude of the fluctuations
            List<Type> fluctuationMagnitude_;

        //- Eckert number used if in Eckert fluctuation magnitude mode.
            List<scalar> EckertNumber_;

        //- The points of the perturbation cell centers
            List<List<vector> > points_;

        //- The bounding boxes for each perturbation zone
            List<boundBox> zoneBox_;
 
        //- The bounding boxes for each perturbation cell
            List<List<boundBox> > cellBox_;

        //- The indices of the grid cells within a perturbation cell
            List<List<List<label> > > gridCellsInCellBox_;

        //- The perturbation cell fluctuation field
            List<List<Type> > fluctuations_;

        //- The mode that the fluctuations update in time.
        //  Options are:
        //     - "fixedFrequency" Entire perturbation field is updated at the rate dictated by 'updatePeriod'.
        //     - "windAtHeight"   Entire perturbation field is updated together at a rate based on a
        //                        horizontal wind vector sampled at some height and the horizontal 
        //                        size of the perturbation zone (i.e., the time for the flow to traverse
        //                        the perturbation zone at some height.
        //     - "slabLocalWind"  Every horizontal perturbation slab is updated independently at a rate
        //                        based on the horizontal wind vector and horizontal size of the 
        //                        perturbation slab (i.e., local flow through time per slab)
            List<word> updateMode_;

        //- Update period.  If "updateMode" is "fixedFrequency", this is
        //  the update period in seconds.
            List<scalar> updatePeriod_;

        //- This dictates how the field is updated.  Options are:
        //  --direct: the perturbations are directly added to the field.
        //  --sourceTerm: the perturbations are introduced through a source
        //                term in the corresponding governing equations RHS.
            List<word> applicationMode_;

        //- Last time that the fluctuations were updated.
            List<scalar> lastUpdateTime_;

        //- Average flow velocity at some height.
            List<vector> averageVelAtHeight_;
            
        //- Average flow velocity over slabs of perturbation cells.
            List<List<vector> > averageVelOverSlabs_;


    // Private Member Functions

        //- Initialize this class
        void initialize();

        //- Read the sub-dictionaries
        void readSubDict();

        //- Check for necessary inputs.  If they are not there, throw an error.
        void inputChecks();

        //- Check to see if wall distance will be used.  If so, compute it.
        void defineWallDist();

        //- Create the perturbation cells (different than grid cells)
        void createPerturbationCells();

        //- Identify grid cells within a perturbation zone.
        void identifyGridCellsInZone(int m, List<label>& gridCellList);

        //- Identify the grid cells that lie within perturbation cells.
        void identifyGridCellsInCellBox();

        //- Identify grid cells within the perturbation zone that lie within
        //  a certain height zone.
        void identifyGridCellsAtHeight(int m, scalar h, List<label>& gridCellList);

        //- Get average velocity of grid cells identified by identifyGridCellsAtHeight
        void getVelocityAtHeight(int m, scalar h, vector velAvg, vector velMin, vector velMax);

        //- Get the average velocity of horizontal slabs of perturbation cells, so
        //  a list of average velocities with height
        void getVelocityOverSlabs(int m, List<vector>& velAvg, List<vector>& velMin, List<vector>& velMax);

        //- Update the fluctuation values in each perturbation cell.  The second version allows
        //  the update of only horizontal slab at height k.
        void updateCellFluctuations(int m);
        void updateCellFluctuations(int m, int k);

        //- Update the source term such that it will change the field by the
        //  specified fluctuation field in one time step.
        void updateSourceTerm(int m);

        //- Zero the source term in perturbation zone m.
        void zeroSourceTerm(int m);

        //- Update the perturbation field.
        void updatePerturbationField(int m);

        //- Zero the perturbation field.
        void zeroPerturbationField(int m);



public:

    //- Runtime type information
    ClassName("perturbationZone");


    // Constructors
    perturbationZone 
    (
        const IOdictionary& dict,
        volFieldType& field
    );


    // Destructor
    virtual ~perturbationZone();


    // Overall update function that applications can call
    void update();

    // Directly add the perturbations to the field.
    void perturb(volFieldType& fld)
    {
        fld += fieldPerturbation_;
    };

    // Return the source term field.
    volFieldType& source()    
    {
        return source_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "perturbationZone.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
