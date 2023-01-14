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
    addToRunTimeSelectionTable
    (
        functionObject,
        geometricBlendingFactor,
        dictionary
    );
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

    Info << "fldExist = " << fldExist << endl;

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


void Foam::functionObjects::geometricBlendingFactor::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Stabilization blending factor");
    writeCommented(os, "Time");
    writeTabbed(os, "Scheme1");
    writeTabbed(os, "Scheme2");
    writeTabbed(os, "Blended");
    os  << endl;
}


bool Foam::functionObjects::geometricBlendingFactor::calc()
{
    init(false);
    return true;
}


bool Foam::functionObjects::geometricBlendingFactor::init(bool first)
{
  //const auto* residualPtr = mesh_.findObject<IOField<scalar>>(residualName_);

    volScalarField& indicator = this->indicator();

    const volScalarField& gridCellVolumeRef = lookupObject<volScalarField>(gridCellVolumeName_);

    if (gridCellVolume_)
    {
        /*
        if (maxNonOrthogonality_ >= minNonOrthogonality_)
        {
            FatalErrorInFunction
                << "minNonOrthogonality should be larger than "
                << "maxNonOrthogonality."
                << exit(FatalError);
        }
        */

      //indicator = scalar(1);
        forAll(indicator, i)
        {
            indicator[i] = mesh_[i].Cf();
        }
        /*
            max
            (
                indicator,
                min
                (
                    max
                    (
                        scalar(0),
                        (*nonOrthPtr - maxNonOrthogonality_)
                      / (minNonOrthogonality_ - maxNonOrthogonality_)
                    ),
                    scalar(1)
                )
            );
        */

        if (first)
        {
            Log << "    Max grid cell volume :  " << max(gridCellVolumeRef).value()
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

        Log << nl << name() << " execute :" << nl
        << "    scheme 1 cells :  " << nCellsScheme1 << nl
        << "    scheme 2 cells :  " << nCellsScheme2 << nl
        << "    blended cells  :  " << nCellsBlended << nl
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
    writeFile(obr_, name),
    gridCellVolume_(dict.lookupOrDefault<Switch>("useGridCellVolume", false)),
    topoSetSource_(dict.lookupOrDefault<Switch>("useTopoSetSource", false)),

    gridCellVolumeName_
    (
        dict.lookupOrDefault<word>("gridCellVolumeName", "gridCellVolume")
    ),
    topoSetSourceName_
    (
        dict.lookupOrDefault<word>("topoSetSourceName", "topoSetSource")
    ),

    tolerance_(0.001)
{
    Info << "A..." << endl;
   
    read(dict);

    Info << "B..." << endl;

    setResultName(typeName, "");

    Info << "C..." << endl;


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

    Info << "D..." << endl;
    store(resultName_, tfaceBlended);
    Info << "E..." << endl;




    bool gridCellVolumeExist = mesh_.foundObject<volScalarField>(gridCellVolumeName_);
    Info << "gridCellVolumeExist = " << gridCellVolumeExist << endl;

    if (!gridCellVolumeExist)
    {
        volScalarField* gridCellVolume
        (
            new volScalarField
            (
                IOobject
                (
                    gridCellVolumeName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
              //mesh_.V(),
                dimensionedScalar("volume",dimVol,Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        mesh_.objectRegistry::store(gridCellVolume);
    }
    volScalarField& gridCellVolume = mesh_.lookupObjectRef<volScalarField>(gridCellVolumeName_);

  //gridCellVolume = mesh_.V();

    gridCellVolumeExist = mesh_.foundObject<volScalarField>(gridCellVolumeName_);
    Info << "gridCellVolumeExist = " << gridCellVolumeExist << endl;

/*
    const volScalarField& topoSetSourceRef = mesh_.lookupObject<volScalarField>(topoSetSourceName_);

    if (topoSetSource_ && isNull(topoSetSourceRef))
    {
        IOobject fieldHeader
        (
            topoSetSourceName_,
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldHeader.typeHeaderOk<volScalarField>(true))
        {
            volScalarField* vfPtr = new volScalarField(fieldHeader, mesh_);
            mesh_.objectRegistry::store(vfPtr);
        }
        else
        {
            FatalErrorInFunction
                << "Field : " << topoSetSourceName_ << " not found."
                << " The function object will not be used"
                << exit(FatalError);
        }
    }
*/

    if (log)
    {
      //indicator().writeOpt(IOobject::AUTO_WRITE);
//        indicator().writeOpt() = IOobject::AUTO_WRITE;
    }

  //if (writeToFile_)
  //{
  //    writeFileHeader(file());
  //}

    Info << "G..." << endl;
    init(true);
    Info << "H..." << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::geometricBlendingFactor::read
(
    const dictionary& dict
)
{
  //if (fieldExpression::read(dict) && writeFile::read(dict))
    {
        dict.lookupOrDefault<Switch>("useGridCellVolume", gridCellVolume_);
        dict.lookupOrDefault<Switch>("useTopoSetSource", topoSetSource_);

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

        Info<< type() << " " << name() << ":" << nl;
        if (gridCellVolume_)
        {
            Info<< "    Including gridCellVolume between: "
              //<< minNonOrthogonality_ <<  " and " << maxNonOrthogonality_
                << endl;
        }
        if (topoSetSource_)
        {
            Info<< "    Including topoSetSource between: "
              //<< minGradCc_ << " and " << maxGradCc_
                << endl;
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

/*
    if (writeToFile_)
    {
        writeCurrentTime(file());

        file()
            << tab << nCellsScheme1
            << tab << nCellsScheme2
            << tab << nCellsBlended
            << endl;
    }
*/

    return true;
}


// ************************************************************************* //
