/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceClean

Group
    grpSurfaceUtilities

Description
    Utility to clean surfaces.

    Current functionality
    - removes baffles
    - collapses small edges, removing triangles.
    - converts sliver triangles into split edges by projecting point onto
      base of triangle.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"

#include "collapseBase.H"
#include "collapseEdge.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Clean surface by removing baffles, sliver faces,"
        " collapsing small edges, etc."
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("length", "The min length");
    argList::addArgument("quality", "The min quality");
    argList::addArgument("output", "The output surface file");

    argList::addBoolOption
    (
        "noClean",
        "Suppress surface checking/cleanup on the input surface"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Input geometry scaling factor"
    );
    argList args(argc, argv);

    const fileName inFileName = args[1];
    const scalar minLen = args.get<scalar>(2);
    const scalar minQuality = args.get<scalar>(3);
    const fileName outFileName = args[4];

    Info<< "Reading surface " << inFileName << nl
        << "Collapsing all triangles with" << nl
        << "    edges or heights < " << minLen << nl
        << "    quality          < " << minQuality << nl
        << "Writing result to " << outFileName << nl << endl;


    Info<< "Reading surface from " << inFileName << " ..." << nl << endl;

    triSurface surf
    (
        inFileName,
        args.opt<scalar>("scale", -1)
    );
    surf.writeStats(Info);

    if (!args.found("noClean"))
    {
        Info<< "Removing duplicate and illegal triangles ..." << nl << endl;
        surf.cleanup(true);
    }

    Info<< "Collapsing triangles to edges ..." << nl << endl;

    while (true)
    {
        const label nEdgeCollapse = collapseEdge(surf, minLen);

        if (nEdgeCollapse == 0)
        {
            break;
        }
    }
    while (true)
    {
        const label nSplitEdge = collapseBase(surf, minLen, minQuality);

        if (nSplitEdge == 0)
        {
            break;
        }
    }

    Info<< nl
        << "Resulting surface:" << endl;
    surf.writeStats(Info);

    Info<< nl
        << "Writing refined surface to " << outFileName << " ..." << endl;
    surf.write(outFileName);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
