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

#include "spongeLayer.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spongeLayer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::spongeLayer::readSingleSpongeSubdict_(int s)
{
    // Get current sponge and read the associated subdict
    currentSponge_ = spongesList_[s];
    dictionary currentSpongeDict(spongeDict_.subDict(currentSponge_));

    // Type of layer 
    type_ = currentSpongeDict.lookupOrDefault<word>("type","none");

    // Sponge layer width
    if ( currentSpongeDict.found("width") )
    {
        // Sponge layer has constant width
        width_ = currentSpongeDict.lookupOrDefault<scalar>("width",5000.0);
    }
    else
    {
        // Sponge layer has variable width. Read table and interpolate
        List<List<scalar>> widthTable( currentSpongeDict.lookup("widthTable") );
        forAll(widthTable, i)
        {
            widthTableTime_.append( widthTable[i][0] );
            widthTableWidth_.append( widthTable[i][1] );
        }

        width_ = interpolateXY(runTime_.value(),widthTableTime_, widthTableWidth_);

        needsLocationUpdating_ = true;
    }

    // Sponge layer parameters
    if ( currentSpongeDict.found("patch") )
    {
        // Standard at-domain-boundary layer

        // Patch to apply the layer
        patch_ = currentSpongeDict.lookupOrDefault<word>("patch","null");

        // Set standard parameters
        const boundBox& bb = mesh_.bounds();
        if (patch_ == "upper")
        {
            startLocation_ = bb.max()[2] - width_;
            coordIndex_ = 2;
            direction_ = "stepUp";
        }
        else if (patch_ == "north")
        {
            startLocation_ = bb.max()[1] - width_;
            coordIndex_ = 1;
            direction_ = "stepUp";
        }
        else if (patch_ == "south")
        {
            startLocation_ = bb.min()[0];
            coordIndex_ = 1;
            direction_ = "stepDown";
        }
        else if (patch_ == "west")
        {
            startLocation_ = bb.min()[0];
            coordIndex_ = 0;
            direction_ = "stepDown";
        }
        else if (patch_ == "east")
        {
            startLocation_ = bb.max()[0] - width_;
            coordIndex_ = 0;
            direction_ = "stepUp";
        }
        else
            Info << "ERROR: patch does not exist or it is not valid." << endl;

    }
    else
    {
        // Non-standard (generic) layer specification

        // Sponge layer start location
        startLocation_ = currentSpongeDict.lookupOrDefault<scalar>("startLocation",0.0);

        // Coordinate index
        coordIndex_ = currentSpongeDict.lookupOrDefault<label>("coordIndex",2);

        // Step up or step down
        direction_ = currentSpongeDict.lookupOrDefault<word>("direction","stepUp");
    }

    // Maximum viscosity
    viscosityMax_ = currentSpongeDict.lookupOrDefault<scalar>("viscosityMax",0.01); 
    
    // Components to actively damp
    dampingComp_ = currentSpongeDict.lookupOrDefault<word>("dampingComp","vertical");

    // Create sponge layer reference velocity
    Ux_ = currentSpongeDict.lookupOrDefault<scalar>("Ux",0.0);
    Uy_ = currentSpongeDict.lookupOrDefault<scalar>("Uy",0.0);

}

void Foam::spongeLayer::update()
{
    // Compute the sponge layer damping force, solving an equation for each type of layer 
    // This approach is contrasted by solving an equation for each layer; however, a viscosity
    // field already summed all the layers of that type

    // First, update the widths of the layers by updating the viscosity/tau fields
    if ( needsLocationUpdating_ )
    {
        updateSpongeLocation_();
    }

    // Nest, calculate the damping force based on the new viscosity/tau fields
    if ( spongesList_.size() == 0 )
        return;

    vector Uref;
    Uref.x() = Ux_;
    Uref.y() = Uy_;
    Uref.z() = 0.0;
    Uref_ = dimensionedVector("Uref", dimensionSet(0, 1, -1, 0, 0, 0, 0), Uref);

    volVectorField Unew= 1*U_;
    forAll(Unew,cellI)
    {
        Unew[cellI].x()=0;
        Unew[cellI].y()=0;
    }

    // "Rayleigh" and "vertical"
    bodyForce_ = - tauV_ * Unew;

    // "Rayleigh" and "horizontal"
    bodyForce_ += tauH_ * (Uref_ - U_);

    // "viscous" and "vertical"
    bodyForce_ += fvc::laplacian(viscosityV_,Unew);

    // "viscous" and "horizontal"
    bodyForce_ += fvc::laplacian(viscosityH_,U_);

}

void Foam::spongeLayer::updateSpongeLocation_()
{
    clearViscosityFields_();

    // Update the viscosity/tau fields
    forAll(spongesList_, s)
    {
        // Read the subdict of the current sponge
        readSingleSpongeSubdict_(s);

        // Add to the appropriate viscosity/tau field
        addSponge_(s);
    }

    adjustOverlappingVisc_();
}


void Foam::spongeLayer::clearViscosityFields_()
{
    tauV_ = 0.0 * tauV_;
    tauH_ = 0.0 * tauH_;
    viscosityV_ = 0.0 * viscosityV_;
    viscosityH_ = 0.0 * viscosityH_;
}


void Foam::spongeLayer::addSponge_(int s)
{
    // Print information
    if ( runTime_.timeIndex() == 0 ) // before time loop
    {
        char xyzdir='?';
        switch(coordIndex_){
            case 0: xyzdir='x'; break;
            case 1: xyzdir='y'; break;
            case 2: xyzdir='z'; break;
        }

        Info << "Adding " << currentSponge_  << ": " << type_ << " damping between ";
        Info << startLocation_ << " and " << startLocation_+width_ << " (" << direction_;
        Info << ") in " << xyzdir << " direction; " << dampingComp_ << " damping with ";
        Info << viscosityMax_ << " max damping coefficient" << endl;
    }
    else
    {
        Info << "Sponge Layers: Updating " << currentSponge_ << " damping layer location to ";
        Info << startLocation_ << " to " << startLocation_+width_ << "." << endl;
    }

    // Set viscosity to cosine profile between startLocation and startLocation+width,
    // For step up:   zero below startLocation and one  above startLocation+width
    // For step down: one  below startLocation and zero above startLocation+width
    // viscosity is the percentage of the viscosityMax set in the dictionary
    scalar fact = 1.0; //stepUp
    if (direction_ == "stepDown")
    {
        fact = -1.0;
    }

    // Get the current layer's type viscosity field from previously added layers
    if (type_ == "Rayleigh" && dampingComp_ == "vertical")
    {
        viscosity_.dimensions().reset(dimensionSet(0, 0, -1, 0, 0, 0, 0));
        viscosity_ = tauV_;
    }
    else if (type_ == "Rayleigh" && dampingComp_ == "horizontal")
    {
        viscosity_.dimensions().reset(dimensionSet(0, 0, -1, 0, 0, 0, 0));
        viscosity_ = tauH_;
    }
    else if (type_ == "viscous"  && dampingComp_ == "vertical")
    {
        viscosity_.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));
        viscosity_ = viscosityV_;
    }
    else if (type_ == "viscous"  && dampingComp_ == "horizontal")
    {
        viscosity_.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));
        viscosity_ = viscosityH_;
    }
    else
    {
        Info << "^ The layer above has not been adeed. Invalid or \"none\" type added." << endl;
    }

    scalar temp;
    forAll(mesh_.cells(),cellI)
    {
        scalar loc = mesh_.C()[cellI][coordIndex_];

        temp  = (loc<=startLocation_) * (1.0 - fact);
        temp += ((loc>startLocation_) && (loc<startLocation_+width_)) *
            (
             1.0 - fact * Foam::cos
             (
              Foam::constant::mathematical::pi * (loc - startLocation_)/width_
             )
            );
        temp += (loc>=startLocation_+width_) * (1.0 + fact);
        temp *= 0.5 * viscosityMax_;
        if (temp > viscosity_[cellI])
            viscosity_[cellI] = temp;
    }

    forAll(viscosity_.boundaryField(),i)
    {
        if ( !mesh_.boundary()[i].coupled() )
        {
            forAll(viscosity_.boundaryField()[i],j)
            {
                scalar loc = mesh_.boundary()[i].Cf()[j][coordIndex_];
                temp  = (loc<=startLocation_) * (1.0 - fact);
                temp += ((loc>startLocation_) && (loc<startLocation_+width_)) *
                    (
                     1.0 - fact * Foam::cos
                     (
                      Foam::constant::mathematical::pi * (loc - startLocation_)/width_
                     )

                    );
                temp += (loc>=startLocation_+width_) * (1.0 + fact);
                temp *= 0.5 * viscosityMax_;
                if (temp > viscosity_.boundaryField()[i][j])
                    viscosity_.boundaryFieldRef()[i][j] = temp;
            }
        }
    }

    // Add the viscosity into the correct variable
    if      (type_ == "Rayleigh" && dampingComp_ == "vertical")
    {
            tauV_ = viscosity_;
    }
    else if (type_ == "Rayleigh" && dampingComp_ == "horizontal")
    {
            tauH_ = viscosity_;
    }
    else if (type_ == "viscous"  && dampingComp_ == "vertical" )
    {
            viscosityV_ = viscosity_;
    }
    else if (type_ == "viscous"  && dampingComp_ == "horizontal")
    {
            viscosityH_ = viscosity_;
    }
    
}


void Foam::spongeLayer::adjustOverlappingVisc_()
{
    // When two viscosity fields (vicosity*_ or tau*_) overlap, only the highest one
    // should be considered for the bodyForce computation. A linear sum will result 
    // in a overestimation of the force field

    scalar maxVisc;
    forAll(mesh_.cells(),cellI)
    {
        maxVisc = max(max(viscosityH_[cellI], viscosityV_[cellI]),max(tauH_[cellI], tauV_[cellI]));
        //gMax for parallel?

        if ( maxVisc == tauV_[cellI] )
        {
            viscosityH_[cellI]=0.0; viscosityV_[cellI]=0.0; tauH_[cellI]=0.0;
        }
        else if ( maxVisc == tauH_[cellI] )
        {
            viscosityH_[cellI]=0.0; viscosityV_[cellI]=0.0; tauV_[cellI]=0.0;
        }
        else if ( maxVisc == viscosityV_[cellI] )
        {  
            viscosityH_[cellI]=0.0; tauH_[cellI]=0.0; tauV_[cellI]=0.0;
        }
        else if ( maxVisc == viscosityH_[cellI] )
        {
            viscosityV_[cellI]=0.0; tauH_[cellI]=0.0; tauV_[cellI]=0.0;
        }
    }

    forAll(viscosity_.boundaryField(),i)
    {
        if ( !mesh_.boundary()[i].coupled() )
        {
            forAll(viscosity_.boundaryField()[i],j)
            {
                maxVisc = max(max(viscosityH_.boundaryField()[i][j], viscosityV_.boundaryField()[i][j]),
                              max(tauH_.boundaryField()[i][j], tauV_.boundaryField()[i][j])   );
                
                if ( maxVisc == tauV_.boundaryField()[i][j] )
                {
                    viscosityH_.boundaryFieldRef()[i][j]=0.0;
                    viscosityV_.boundaryFieldRef()[i][j]=0.0;
                    tauH_.boundaryFieldRef()[i][j]=0.0;
                }
                else if ( maxVisc == tauH_.boundaryField()[i][j] )
                {
                    viscosityH_.boundaryFieldRef()[i][j]=0.0;
                    viscosityV_.boundaryFieldRef()[i][j]=0.0;
                    tauV_.boundaryFieldRef()[i][j]=0.0;
                }
                else if ( maxVisc == viscosityV_.boundaryField()[i][j] )
                {  
                    viscosityH_.boundaryFieldRef()[i][j]=0.0;
                    tauH_.boundaryFieldRef()[i][j]=0.0;
                    tauV_.boundaryFieldRef()[i][j]=0.0;
                }
                else if ( maxVisc == viscosityH_.boundaryField()[i][j] )
                {
                    viscosityV_.boundaryFieldRef()[i][j]=0.0;
                    tauH_.boundaryFieldRef()[i][j]=0.0;
                    tauV_.boundaryFieldRef()[i][j]=0.0;
                }
            }
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spongeLayer::spongeLayer
(
    const word& name,
    const volVectorField& U
)
:
    // Set name
    name_(name),

    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity field
    U_(U),

    // Initialize the reference velocity field
    Uref_
    (
        IOobject
        (
            name_ & "Uref",
            runTime_.constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedVector("Uref_", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
    ),

    // Initialize the viscosity fields
    // `viscosity` if viscous, `tau` if Rayleigh
    viscosity_
    (
        IOobject
        (
            name_ & "viscosity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE //NO_WRITE
        ),
        mesh_,
        dimensionedScalar("viscosity_", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    viscosityH_
    (
        IOobject
        (
            name_ & "viscosityH",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE //NO_WRITE
        ),
        mesh_,
        dimensionedScalar("viscosityH_", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    viscosityV_
    (
        IOobject
        (
            name_ & "viscosityV",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE //NO_WRITE
        ),
        mesh_,
        dimensionedScalar("viscosityV_", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    tauH_
    (
        IOobject
        (
            name_ & "tauH",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh_,
        dimensionedScalar("tauH_", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    tauV_
    (
        IOobject
        (
            name_ & "tauV",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh_,
        dimensionedScalar("tauV_", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),

    // Initialize the body force field
    bodyForce_
    (
        IOobject
        (
            name_ & "Force",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
    )


{
    // Define dictionary with input data
    IOdictionary ABLProperties
    (
        IOobject
        (
            "ABLProperties",
            runTime_.time().constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    spongeDict_ = ABLProperties.subOrEmptyDict(name_);
    spongesList_ = spongeDict_.toc();

    if ( spongesList_.size() == 0 )
    {
        Info << "No sponge layers specified in ABLProperties. Skipping." << endl;
    }
  
    // Unless there is a sponge layer with time-varying width, their location will
    // not be updated during the main time loop. This variable will be set to true
    // if a widthTable is encontered and the loop on all sponges will be performed
    // at every time step.
    needsLocationUpdating_ = false;
        
    // Loop over all sponge subdictionaries
    forAll(spongesList_, s)
    {
        // Read the subdict of the current sponge
        readSingleSpongeSubdict_(s);

        // Add to the appropriate viscosity/tau field
        addSponge_(s);
    }

    adjustOverlappingVisc_();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spongeLayer::~spongeLayer()
{}


// ************************************************************************* //
