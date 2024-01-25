/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "geometricBlendingFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(geometricBlendingFactor, 0);
    addToRunTimeSelectionTable(functionObject, geometricBlendingFactor, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Gives the pointer to the indicator field (i.e., the field that indicates
// what the blending function should be).
Foam::volScalarField&
Foam::functionObjects::geometricBlendingFactor::indicator()
{
    // The indicator field name is "blendedIndicator" plus the
    // name of the field blending is applied to.  For example
    // "blendedIndicatorU" for the velocity field.
    const word fldName = "blendedIndicator" + fieldName_;

    // Check to see if the indicator field is already in the
    // registry.
    bool fldExist = mesh_.foundObject<volScalarField>(fldName);

    // Create a reference to the indicator field.  If the field is not in
    // the registry, a null object will be returned.
    if (!fldExist)
    {
        volScalarField* fld
        (
            new volScalarField
            (
                IOobject
                (
                    "blendedIndicator" + fieldName_,
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero",dimless,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(fld);
    }
    volScalarField& fld = mesh_.lookupObjectRef<volScalarField>(fldName);


    return fld;
}


void Foam::functionObjects::geometricBlendingFactor::calcStats
(
    label& nCellsScheme1,
    label& nCellsScheme2,
    label& nCellsBlended
) const
{
    const word indicatorName = "blendedIndicator" + fieldName_;
    const volScalarField& indicator = mesh_.lookupObject<volScalarField>(indicatorName);
    if (isNull(indicator))
    {
        FatalErrorInFunction
            << "Indicator field not set"
            << abort(FatalError);
    }

    for (const scalar i : indicator)
    {
        if (i < tolerance_)
        {
            ++nCellsScheme2;
        }
        else if (i > (1 - tolerance_))
        {
            ++nCellsScheme1;
        }
        else
        {
            ++nCellsBlended;
        }
    }

    reduce(nCellsScheme1, sumOp<label>());
    reduce(nCellsScheme2, sumOp<label>());
    reduce(nCellsBlended, sumOp<label>());
}


void Foam::functionObjects::geometricBlendingFactor::writeFileHeader(const label i)
{
    writeHeader(file(), "Stabilization blending factor");
    writeCommented(file(), "Time");
    writeTabbed(file(), "Scheme1");
    writeTabbed(file(), "Scheme2");
    writeTabbed(file(), "Blended");
    file()  << endl;
}


bool Foam::functionObjects::geometricBlendingFactor::calc()
{
    read(dict_);
    init(false);
    return true;
}


bool Foam::functionObjects::geometricBlendingFactor::init(bool first)
{
    volScalarField& indicator = this->indicator();
    indicator = 0.0;

    const volScalarField& gridCellResolutionRef = lookupObject<volScalarField>(gridCellResolutionName_);
    if (gridCellResolution_)
    {
        forAll(indicator, i)
        {
            scalar indicatorCell = gridCellResolutionBlendingTable_->value(gridCellResolutionRef[i]);
            indicator[i] = max(indicator[i],max(min(indicatorCell,scalar(1)),scalar(0)));
        }
        if (first)
        {
            Log << "    Max / Min grid cell resolution:  " << max(gridCellResolutionRef).value() 
                << " / "                                   << min(gridCellResolutionRef).value()
                << endl;
        }
    }

    const volScalarField& gridCellxRef = lookupObject<volScalarField>(gridCellxName_);
    if (gridCellx_)
    {
        forAll(indicator, i)
        {
            scalar indicatorCell = gridCellxBlendingTable_->value(gridCellxRef[i]);
            indicator[i] = max(indicator[i],max(min(indicatorCell,scalar(1)),scalar(0)));
        }
        if (first)
        {
            Log << "    Max / Min grid cell x-value:  " << max(gridCellxRef).value() 
                << " / "                                << min(gridCellxRef).value()
                << endl;
        }
    }

    const volScalarField& gridCellyRef = lookupObject<volScalarField>(gridCellyName_);
    if (gridCelly_)
    {
        forAll(indicator, i)
        {
            scalar indicatorCell = gridCellyBlendingTable_->value(gridCellyRef[i]);
            indicator[i] = max(indicator[i],max(min(indicatorCell,scalar(1)),scalar(0)));
        }
        if (first)
        {
            Log << "    Max / Min grid cell y-value:  " << max(gridCellyRef).value() 
                << " / "                                << min(gridCellyRef).value()
                << endl;
        }
    }

    const volScalarField& gridCellzRef = lookupObject<volScalarField>(gridCellzName_);
    if (gridCellz_)
    {
        forAll(indicator, i)
        {
            scalar indicatorCell = gridCellzBlendingTable_->value(gridCellzRef[i]);
            indicator[i] = max(indicator[i],max(min(indicatorCell,scalar(1)),scalar(0)));
        }
        if (first)
        {
            Log << "    Max / Min grid cell z-value:  " << max(gridCellzRef).value() 
                << " / "                                << min(gridCellzRef).value()
                << endl;
        }
    }


    if (first)
    {
        Log << nl;
    }

    if (log)
    {
        label nCellsScheme1 = 0;
        label nCellsScheme2 = 0;
        label nCellsBlended = 0;

        calcStats(nCellsScheme1, nCellsScheme2, nCellsBlended);

        Log << nl << name() << " execute:" << nl
        << "    scheme 1 cells:  " << nCellsScheme1 << nl
        << "    scheme 2 cells:  " << nCellsScheme2 << nl
        << "    blended cells:   " << nCellsBlended << nl
        << endl;
    }

    indicator.correctBoundaryConditions();
    indicator.min(1.0);
    indicator.max(0.0);

    // Update the blended surface field
    auto& surBlended = mesh_.lookupObjectRef<surfaceScalarField>(resultName_);

    surBlended = fvc::interpolate(indicator);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::geometricBlendingFactor::geometricBlendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    logFiles(obr_, name),

    gridCellResolution_(dict.lookupOrDefault<Switch>("useGridCellResolution", false)),
    gridCellx_(dict.lookupOrDefault<Switch>("useGridCellxValue", false)),
    gridCelly_(dict.lookupOrDefault<Switch>("useGridCellyValue", false)),
    gridCellz_(dict.lookupOrDefault<Switch>("useGridCellzValue", false)),

    gridCellResolutionBlendingTable_(),
    gridCellxBlendingTable_(),
    gridCellyBlendingTable_(),
    gridCellzBlendingTable_(),

    gridCellResolutionName_(dict.lookupOrDefault<word>("gridCellResolutionName", "gridCellResolution")),
    gridCellxName_(dict.lookupOrDefault<word>("gridCellxValueName", "gridCellxValue")),
    gridCellyName_(dict.lookupOrDefault<word>("gridCellyValueName", "gridCellyValue")),
    gridCellzName_(dict.lookupOrDefault<word>("gridCellzValueName", "gridCellzValue")),

    tolerance_(0.001),

    dict_(dict)
{

    // Read the dictionary input for this function object.
    read(dict_);

    // Set the name of the blending factor field based.
    setResultName(typeName, "");

    // Set the name of the summary file written in postProcessing.
    resetName(typeName);
    createFiles();

    // Create an IO object that is the blending field, which is a grid cell
    // face field (because the div scheme interpolates to faces) and store
    // it in the registry to be able to access later.
    tmp<surfaceScalarField> tfaceBlended
    (
        new surfaceScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, Zero)
        )
    );
    store(resultName_, tfaceBlended);


    // Check to see if the grid cell volume field exists.  If it does not, create
    // it and put it in the IO object registry.  Then set the gridCellVolume variable
    // as a reference to this field.
    bool gridCellResolutionExist = mesh_.foundObject<volScalarField>(gridCellResolutionName_);

    if (!gridCellResolutionExist)
    {
        volScalarField* gridCellResolution
        (
            new volScalarField
            (
                IOobject
                (
                    gridCellResolutionName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("resolution",dimLength,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(gridCellResolution);
    }
    volScalarField& gridCellResolution = mesh_.lookupObjectRef<volScalarField>(gridCellResolutionName_);

    // Set the gridCellResolution field to the actual grid cell resolution in case it has not already
    // been set, which would the the usual case.
    forAll(mesh_.V(), i)
    {
        gridCellResolution[i] = cbrt(mesh_.V()[i]);
    }
    gridCellResolution.correctBoundaryConditions();



    // Check to see if the grid cell x-value field exists.  If it does not, create
    // it and put it in the IO object registry.  Then set the gridCellz variable
    // as a reference to this field.
    bool gridCellxExist = mesh_.foundObject<volScalarField>(gridCellxName_);

    if (!gridCellxExist)
    {
        volScalarField* gridCellx
        (
            new volScalarField
            (
                IOobject
                (
                    gridCellxName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("xValue",dimLength,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(gridCellx);
    }
    volScalarField& gridCellx = mesh_.lookupObjectRef<volScalarField>(gridCellxName_);

    // Set the gridCellz field to the actual grid cell z-value in case it has not already
    // been set, which would the the usual case.
    forAll(mesh_.C(), i)
    {
        gridCellx[i] = mesh_.C()[i].x();
    }
    gridCellx.correctBoundaryConditions();


    // Check to see if the grid cell y-value field exists.  If it does not, create
    // it and put it in the IO object registry.  Then set the gridCellz variable
    // as a reference to this field.
    bool gridCellyExist = mesh_.foundObject<volScalarField>(gridCellyName_);

    if (!gridCellyExist)
    {
        volScalarField* gridCelly
        (
            new volScalarField
            (
                IOobject
                (
                    gridCellyName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("yValue",dimLength,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(gridCelly);
    }
    volScalarField& gridCelly = mesh_.lookupObjectRef<volScalarField>(gridCellyName_);

    // Set the gridCellz field to the actual grid cell z-value in case it has not already
    // been set, which would the the usual case.
    forAll(mesh_.C(), i)
    {
        gridCelly[i] = mesh_.C()[i].y();
    }
    gridCelly.correctBoundaryConditions();


    // Check to see if the grid cell z-value field exists.  If it does not, create
    // it and put it in the IO object registry.  Then set the gridCellz variable
    // as a reference to this field.
    bool gridCellzExist = mesh_.foundObject<volScalarField>(gridCellzName_);

    if (!gridCellzExist)
    {
        volScalarField* gridCellz
        (
            new volScalarField
            (
                IOobject
                (
                    gridCellzName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zValue",dimLength,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(gridCellz);
    }
    volScalarField& gridCellz = mesh_.lookupObjectRef<volScalarField>(gridCellzName_);

    // Set the gridCellz field to the actual grid cell z-value in case it has not already
    // been set, which would the the usual case.
    forAll(mesh_.C(), i)
    {
        gridCellz[i] = mesh_.C()[i].z();
    }
    gridCellz.correctBoundaryConditions();



    // If the log option is true, then set the indicator field to write at 
    // write time.
    if (log)
    {
        indicator().writeOpt() = IOobject::AUTO_WRITE;
    }

    // Call the initialization function for the first time.
    init(true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::geometricBlendingFactor::read
(
    const dictionary& dict
)
{
    if (fieldExpression::read(dict))
    {
        gridCellResolution_ = dict.lookupOrDefault<Switch>("useGridCellResolution", false);
        gridCellx_ = dict.lookupOrDefault<Switch>("useGridCellxValue", false);
        gridCelly_ = dict.lookupOrDefault<Switch>("useGridCellyValue", false);
        gridCellz_ = dict.lookupOrDefault<Switch>("useGridCellzValue", false);

        gridCellResolutionName_ = dict.lookupOrDefault<word>("gridCellResolutionName", "gridCellResolution");
        gridCellxName_ = dict.lookupOrDefault<word>("gridCellxValueName", "gridCellxValue");
        gridCellyName_ = dict.lookupOrDefault<word>("gridCellyValueName", "gridCellyValue");
        gridCellzName_ = dict.lookupOrDefault<word>("gridCellzValueName", "gridCellzValue");

        if (gridCellResolution_)
        {
            gridCellResolutionBlendingTable_ = Function1<scalar>::New("gridCellResolutionBlendingTable",dict);
        }
        if (gridCellx_)
        {
            gridCellxBlendingTable_ = Function1<scalar>::New("gridCellxValueBlendingTable",dict);
        }
        if (gridCelly_)
        {
            gridCellyBlendingTable_ = Function1<scalar>::New("gridCellyValueBlendingTable",dict);
        }
        if (gridCellz_)
        {
            gridCellzBlendingTable_ = Function1<scalar>::New("gridCellzValueBlendingTable",dict);
        }
    
        tolerance_ = 0.001;
        if
        (
            dict.readIfPresent("tolerance", tolerance_)
         && (tolerance_ < 0 || tolerance_ > 1)
        )
        {
            FatalErrorInFunction
                << "tolerance must be in the range 0 to 1.  Supplied value: "
                << tolerance_ << exit(FatalError);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::geometricBlendingFactor::write()
{
    label nCellsScheme1 = 0;
    label nCellsScheme2 = 0;
    label nCellsBlended = 0;

    calcStats(nCellsScheme1, nCellsScheme2, nCellsBlended);

    if (Pstream::master())
    {
        writeTime(file());
        file()
            << tab << nCellsScheme1
            << tab << nCellsScheme2
            << tab << nCellsBlended
            << endl;
    }

    return true;
}


// ************************************************************************* //
