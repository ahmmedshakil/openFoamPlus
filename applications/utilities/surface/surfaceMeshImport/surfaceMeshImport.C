/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2018 OpenCFD Ltd.
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
    surfaceMeshImport

Group
    grpSurfaceUtilities

Description
    Import from various third-party surface formats into surfMesh
    with optional scaling or transformations (rotate/translate)
    on a coordinateSystem.

Usage
    \b surfaceMeshImport inputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface.

      - \par -name \<name\>
        Specify an alternative surface name when writing.

      - \par -scaleIn \<scale\>
        Specify a scaling factor when reading files.

      - \par -scaleOut \<scale\>
        Specify a scaling factor when writing files.

      - \par -dict \<dictionary\>
        Use alternative dictionary for constant/coordinateSystems.

      - \par -from \<coordinateSystem\>
        Specify a coordinate system when reading files.

      - \par -to \<coordinateSystem\>
        Specify a coordinate system when writing files.

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "coordinateSystems.H"
#include "cartesianCS.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Import from various third-party surface formats into surfMesh"
    );

    argList::noParallel();
    argList::addArgument("surface", "The input surface file");

    argList::addBoolOption
    (
        "clean",
        "Perform some surface checking/cleanup on the input surface"
    );
    argList::addOption
    (
        "name",
        "name",
        "Specify an alternative surface name when writing - "
        "default is 'default'"
    );
    argList::addOption
    (
        "scaleIn",
        "factor",
        "Geometry scaling factor on input - default is 1"
    );
    argList::addOption
    (
        "scaleOut",
        "factor",
        "Geometry scaling factor on output - default is 1"
    );
    argList::addOption("dict", "file", "Use alternative coordinateSystems");

    argList::addOption
    (
        "from",
        "coordinateSystem",
        "Specify a local coordinate system when reading files.",
        true // advanced
    );
    argList::addOption
    (
        "to",
        "coordinateSystem",
        "Specify a local coordinate system when writing files.",
        true // advanced
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // try for the latestTime, but create "constant" as needed
    instantList Times = runTime.times();
    if (Times.size())
    {
        label startTime = Times.size()-1;
        runTime.setTime(Times[startTime], startTime);
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);
    }


    const fileName importName = args[1];
    const word exportName = args.opt<word>("name", "default");

    // check that reading is supported
    if (!MeshedSurface<face>::canRead(importName, true))
    {
        return 1;
    }


    // The coordinate transformations (must be cartesian)
    autoPtr<coordSystem::cartesian> fromCsys;
    autoPtr<coordSystem::cartesian> toCsys;

    if (args.found("from") || args.found("to"))
    {
        IOobject ioCsys = IOobject::selectIO
        (
            IOobject
            (
                coordinateSystems::typeName,
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            args.opt<fileName>("dict", "")
        );

        if (!ioCsys.typeHeaderOk<coordinateSystems>(false))
        {
            FatalErrorInFunction
                << ioCsys.objectPath() << nl
                << exit(FatalError);
        }

        coordinateSystems globalCoords(ioCsys);

        if (args.found("from"))
        {
            const word csName(args["from"]);
            const auto* csPtr = globalCoords.lookupPtr(csName);

            if (!csPtr)
            {
                FatalErrorInFunction
                    << "Cannot find -from " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            fromCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        if (args.found("to"))
        {
            const word csName(args["to"]);
            const auto* csPtr = globalCoords.lookupPtr(csName);

            if (!csPtr)
            {
                FatalErrorInFunction
                    << "Cannot find -to " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            toCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        // Maybe fix this later
        if (fromCsys && toCsys)
        {
            FatalErrorInFunction
                << "Only allowed '-from' or '-to' option at the moment."
                << exit(FatalError);
        }
    }


    MeshedSurface<face> surf(importName);

    if (args.found("clean"))
    {
        surf.cleanup(true);
    }


    scalar scaleIn = 0;
    if (args.readIfPresent("scaleIn", scaleIn) && scaleIn > 0)
    {
        Info<< "scale input " << scaleIn << endl;
        surf.scalePoints(scaleIn);
    }

    if (fromCsys)
    {
        Info<< "move points from coordinate system: "
            << fromCsys->name() << endl;
        tmp<pointField> tpf = fromCsys->localPosition(surf.points());
        surf.movePoints(tpf());
    }

    if (toCsys)
    {
        Info<< "move points to coordinate system: "
            << toCsys->name() << endl;
        tmp<pointField> tpf = toCsys->globalPosition(surf.points());
        surf.movePoints(tpf());
    }

    scalar scaleOut = 0;
    if (args.readIfPresent("scaleOut", scaleOut) && scaleOut > 0)
    {
        Info<< "scale output " << scaleOut << endl;
        surf.scalePoints(scaleOut);
    }

    surfMesh smesh
    (
        IOobject
        (
            exportName,
            runTime.constant(),
            runTime
        ),
        std::move(surf)
    );


    Info<< "writing surfMesh:\n  " << smesh.objectPath() << endl;
    smesh.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
