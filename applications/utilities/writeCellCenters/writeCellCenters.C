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

Application
    writeCellCenters

Description
    Helper utility to output cell center locations for generating complex
    initial fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        fileName outputFile
        (
            mesh.time().path()/"constant"/"cellCenters"
        );
        Info<< "Writing mesh cell centers to " << outputFile << endl;
        OFstream os(outputFile);
        os.format("ascii");
        os << mesh.cellCentres();

        fileName outputCSVFile
        (
            mesh.time().path()/"constant"/"cellCenters.csv"
        );
        Info<< "Writing mesh cell centers to " << outputCSVFile << endl;
        OFstream csv(outputCSVFile);
        csv << "x,y,z" << nl;
        forAll(mesh.cellCentres(), cellI)
        {
            vector cc = mesh.cellCentres()[cellI];
            csv << cc.x() << "," << cc.y() << "," << cc.z() << nl;
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
