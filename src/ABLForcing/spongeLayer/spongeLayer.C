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
        width_ = readTableAndInterpolate_(currentSpongeDict, "widthTable");
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
            startLocation_ = bb.min()[1];
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
    if ( currentSpongeDict.found("dampCoeffMax") )
    {
        // Sponge layer has constant (in time) maximum damping coefficient
        dampCoeffMax_ = currentSpongeDict.lookupOrDefault<scalar>("dampCoeffMax",0.01); 
    }
    else
    {
        // Sponge layer has variable damp coefficient. Read table and interpolate
        dampCoeffMax_ = readTableAndInterpolate_(currentSpongeDict, "dampCoeffMaxTable");
    }

    // Components to actively damp
    dampingComp_ = currentSpongeDict.lookupOrDefault<word>("dampingComp","vertical");

    // Create sponge layer reference velocity
    Ux_ = currentSpongeDict.lookupOrDefault<scalar>("Ux",0.0);
    Uy_ = currentSpongeDict.lookupOrDefault<scalar>("Uy",0.0);

    // Fraction of the layer's width with a cosine distribution
    if ( currentSpongeDict.found("cosFractionTable") )
    {
        // Sponge layer has variable cos fraction. Read table and interpolate
        cosFraction_ = readTableAndInterpolate_(currentSpongeDict, "cosFractionTable");
    }
    else
    {
        // Sponge layer has constant (in time) fraction of cosine over its width
        cosFraction_ = currentSpongeDict.lookupOrDefault<scalar>("cosFraction",1.0);
    }

    // Whether or not vertical filter is enabled
    verticalFilter_ = currentSpongeDict.lookupOrDefault<bool>("vertFilter",false);

    // Bool to determine if absolute height or agl height are used
    useWallDistZ_ = currentSpongeDict.lookupOrDefault<bool>("useWallDist", false);

    // Vertical filter start location
    if ( currentSpongeDict.found("vertFilterStartHeightTable") )
    {
        // Vertical filter has variable start location. Read table and interpolate
        vertFiltStartHeight_ = readTableAndInterpolate_(currentSpongeDict, "vertFilterStartHeightTable");
    }
    else
    {
        // Vertical filter has constant start location
        vertFiltStartHeight_ = currentSpongeDict.lookupOrDefault<scalar>("vertFilterStartHeight",0.0);
    }

    // Vertical filter cosine transition thickness
    if ( currentSpongeDict.found("vertFilterCosThicknessTable") )
    {
        // Vertical filter has variable start location. Read table and interpolate
        vertFiltCosThickness_ = readTableAndInterpolate_(currentSpongeDict, "vertFilterCosThicknessTable");
    }
    else
    {
        // Vertical filter has constant start location
        vertFiltCosThickness_ = currentSpongeDict.lookupOrDefault<scalar>("vertFilterCosThickness",100.0);
    }

}

scalar Foam::spongeLayer::readTableAndInterpolate_(dictionary& dict, word parameterTableName)
{
    scalarField paramTableTime_;
    scalarField paramTableParam_;

    List<List<scalar>> paramTable( dict.lookup(parameterTableName) );
    forAll(paramTable, i)
    {
        paramTableTime_.append(  paramTable[i][0] );
        paramTableParam_.append( paramTable[i][1] );
    }

    needsUpdating_ = true;

    return interpolateXY(runTime_.value(),paramTableTime_, paramTableParam_);
}

void Foam::spongeLayer::update()
{
    // Compute the sponge layer damping force, solving an equation for each type of layer 
    // This approach is contrasted by solving an equation for each layer; however, a viscosity
    // field already summed all the layers of that type

    // First, update the parameters of the layers by updating the viscosity/tau fields
    if ( needsUpdating_ )
    {
        // Clear old viscosity field and get updated one
        clearViscosityFields_();
        getViscosityField_();
    }

    // Next, calculate the damping force based on the new viscosity/tau fields
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

void Foam::spongeLayer::getViscosityField_()
{
    // Main loop to read and apply all individual layer parameters

    // Loop over all sponges subdictionaries and get viscosity field
    forAll(spongesList_, s)
    {
        // Read the subdict of the current sponge
        readSingleSpongeSubdict_(s);

        // Determine current sponge's viscosity
        calculateCurrentSpongeViscosity_();

        // Add to the appropriate viscosity/tau field
        addSponge_();
    }

    // Adjust overlap between viscosity{H.V} and tau{H,V}
    adjustOverlappingVisc_();
}


void Foam::spongeLayer::clearViscosityFields_()
{
    tauV_ = 0.0 * tauV_;
    tauH_ = 0.0 * tauH_;
    viscosityV_ = 0.0 * viscosityV_;
    viscosityH_ = 0.0 * viscosityH_;
}

void Foam::spongeLayer::calculateCurrentSpongeViscosity_()
{
    scalar temp;
    scalar start, widthcos, endcos;

    // Set viscosity to cosine profile between startLocation and startLocation+width,
    // For step up:   zero below startLocation and one  above startLocation+width
    // For step down: one  below startLocation and zero above startLocation+width
    // viscosity is the percentage of the dampCoeffMax set in the dictionary
    // Some part of the layer (indicated by cosFraction) will have a constant
    // maximum damping, that is, the cosine will only take place on cosFraction
    // portion of the layer
    scalar fact = 1.0; //stepUp
    if (direction_ == "stepDown")
    {
        fact = -1.0;
    }


    start = startLocation_ + (1.0-fact)/2 * width_*(1-cosFraction_);
    widthcos = max(width_*cosFraction_, SMALL);
    endcos = start+widthcos;

    forAll(mesh_.cells(),cellI)
    {
        scalar loc = mesh_.C()[cellI][coordIndex_];

        temp  = (loc<=start) * (1.0 - fact);
        temp += ((loc>start) && (loc<endcos)) *
            (
             1.0 - fact * Foam::cos( Foam::constant::mathematical::pi * (loc - start)/widthcos )
            );
        temp += (loc>=endcos) * (1.0 + fact);
        temp *= 0.5 * dampCoeffMax_;
        currentViscosity_[cellI] = temp;
    }

    forAll(currentViscosity_.boundaryField(),i)
    {
        if ( !mesh_.boundary()[i].coupled() )
        {
            forAll(currentViscosity_.boundaryField()[i],j)
            {
                scalar loc = mesh_.boundary()[i].Cf()[j][coordIndex_];
                temp  = (loc<=start) * (1.0 - fact);
                temp += ((loc>start) && (loc<endcos)) *
                    (
                     1.0 - fact * Foam::cos ( Foam::constant::mathematical::pi * (loc - start)/widthcos )
                    );
                temp += (loc>=endcos) * (1.0 + fact);
                temp *= 0.5 * dampCoeffMax_;
                currentViscosity_.boundaryFieldRef()[i][j] = temp;
            }
        }
    }

    // Apply the vertical filter if needed
    if (verticalFilter_)
    {
        if ( coordIndex_ == 2)
        {
            Info << "    WARNING: horizontal layers are incompatible with vertical filter" << endl;
        }
        else
        {
            applyVerticalFilter_();
        }
    }


}

void Foam::spongeLayer::applyVerticalFilter_()
{
    // Compute wall distance
    wallDist d(mesh_);
    surfaceScalarField dFace = fvc::interpolate(d.y());

    scalar height, temp;
    scalar start, widthcos, endcos;

    start = vertFiltStartHeight_;
    widthcos = max(vertFiltCosThickness_, SMALL);
    endcos = start+widthcos;

    forAll(mesh_.cells(),cellI)
    {
        if ( useWallDistZ_ )
            height = d.y()[cellI];
        else
            height = mesh_.C()[cellI][2];

        temp  = (height<=start) * 0;
        temp += ((height>start) && (height<endcos)) * 0.5 *
            (
             1.0 - 1.0 * Foam::cos( Foam::constant::mathematical::pi * (height - start)/widthcos )
            );
        temp += (height>=endcos);
        if (currentViscosity_[cellI] > 0)
            currentViscosity_[cellI] *= temp;
    }

    forAll(currentViscosity_.boundaryField(),i)
    {
        if ( !mesh_.boundary()[i].coupled() )
        {
            forAll(currentViscosity_.boundaryField()[i],j)
            {
                if ( useWallDistZ_ )
                    height = dFace.boundaryField()[i][j];
                else
                    height = mesh_.boundary()[i].Cf()[j][2];

                temp  = (height<=start) * 0;
                temp += ((height>start) && (height<endcos)) * 0.5 *
                    (
                     1.0 - 1.0 * Foam::cos ( Foam::constant::mathematical::pi * (height - start)/widthcos )
                    );
                temp += (height>=endcos);
                if (currentViscosity_.boundaryField()[i][j] > 0)
                    currentViscosity_.boundaryFieldRef()[i][j] *= temp;
            }
        }
    }

}


void Foam::spongeLayer::addSponge_()
{
    // Adds the current sponge's viscosity to the appropriate field (viscosity{H,V}, or tau{H,V})

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
        Info << dampCoeffMax_ << " max damping coefficient" << endl;
    }
    else
    {
        Info << "Sponge Layers: Updating " << currentSponge_ << " damping layer location to ";
        Info << startLocation_ << " to " << startLocation_+width_ << " m";
        if ( cosFraction_ == 1.0 )
            Info << "." << endl;
        else
            Info<< "; outer " << (1-cosFraction_)*width_ << " m of constant maximum damping coefficient (" << cosFraction_*100 << "%)"  <<endl;

        Info << "               Updating " << currentSponge_ << " damping layer maximum damping coefficient to ";
        Info << dampCoeffMax_ << endl;

        if ( verticalFilter_ ){
            Info << "               Updating " << currentSponge_ << " vertical filter. No damping applied below ";
            Info << vertFiltStartHeight_ << "; regular damping above " << vertFiltStartHeight_+vertFiltCosThickness_;
            if ( useWallDistZ_ )
                Info << " (agl height)" << endl;
            else
                Info << " (absolute height)" << endl;
        }
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


    // Now, add the current sponge's viscosity, mindful of overlaps in internalField
    // The boundaries are usually not needed, but leaving here in the rare case of a custom
    // layer that is in the middle of the domain but touches some boundaries
    scalar maxVisc;
    forAll (mesh_.cells(), cellI)
    {
        maxVisc = max(viscosity_[cellI], currentViscosity_[cellI]);
        if (maxVisc == currentViscosity_[cellI] )
            viscosity_[cellI] = currentViscosity_[cellI];
    }
    forAll(viscosity_.boundaryField(), i)
    {
        if ( !mesh_.boundary()[i].coupled() )
        {
            forAll(viscosity_.boundaryField()[i],j)
            {
                maxVisc = max(viscosity_.boundaryField()[i][j], currentViscosity_.boundaryField()[i][j]);
                
                if ( maxVisc == currentViscosity_.boundaryField()[i][j] )
                {
                    viscosity_.boundaryFieldRef()[i][j] = currentViscosity_.boundaryField()[i][j];
                }
            }
        }
    }

    
    // Redefine the appropriate viscosity type field with the updated viscosity
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
    // Adjust the overlap between the different types of viscosities (viscosity{H,V},
    // and tau{H,V}. When two of more viscosity fields overlap, only the highest one
    // should be considered for the bodyForce computation. A linear sum will result 
    // in a overestimation of the force field

    scalar maxVisc;
    forAll(mesh_.cells(),cellI)
    {
        maxVisc = max(max(viscosityH_[cellI], viscosityV_[cellI]),max(tauH_[cellI], tauV_[cellI]));

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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("tauV_", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    currentViscosity_
    (
        IOobject
        (
            name_ & "currentViscosity",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("currentViscosity_", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
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
  
    // Unless there is a sponge layer with time-varying parameters, their viscosity
    // field will not be updated during the main time loop. This variable will be 
    // set to true if a time-varying field (identified by the *Table name) is found
    // and thus a new determination of the current overall viscosity field will be
    // performed at every time step.
    needsUpdating_ = false;
        
    // Determine the overall viscosity field where the body force will be applied
    getViscosityField_();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spongeLayer::~spongeLayer()
{}


// ************************************************************************* //
