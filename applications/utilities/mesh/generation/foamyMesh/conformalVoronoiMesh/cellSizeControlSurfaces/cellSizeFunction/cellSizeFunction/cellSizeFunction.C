/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "cellSizeFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellSizeFunction, 0);
    defineRunTimeSelectionTable(cellSizeFunction, dictionary);
}

Foam::scalar Foam::cellSizeFunction::snapToSurfaceTol_ = 1e-10;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeFunction::cellSizeFunction
(
    const word& type,
    const dictionary& cellSizeFunctionDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize,
    const labelList regionIndices
)
:
    dictionary(cellSizeFunctionDict),
    surface_(surface),
    surfaceCellSizeFunction_
    (
        surfaceCellSizeFunction::New
        (
            cellSizeFunctionDict,
            surface,
            defaultCellSize
        )
    ),
    coeffsDict_(optionalSubDict(type + "Coeffs")),
    defaultCellSize_(defaultCellSize),
    regionIndices_(regionIndices),
    sideMode_(),
    priority_
    (
        cellSizeFunctionDict.get<label>("priority", keyType::REGEX_RECURSIVE)
    )
{
    const word mode =
        cellSizeFunctionDict.get<word>("mode", keyType::REGEX_RECURSIVE);

    if (surface_.hasVolumeType())
    {
        if (mode == "inside")
        {
            sideMode_ = smInside;
        }
        else if (mode == "outside")
        {
            sideMode_ = smOutside;
        }
        else if (mode == "bothSides")
        {
            sideMode_ = rmBothsides;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mode, expected: inside, outside or bothSides" << nl
                << exit(FatalError);
        }
    }
    else
    {
        if (mode != "bothSides")
        {
            WarningInFunction
                << "surface does not support volumeType, defaulting mode to "
                << "bothSides."
                << endl;
        }

        sideMode_ = rmBothsides;
    }

    if (debug)
    {
        Info<< nl
            << "Cell size function for surface " << surface.name()
            << ", " << mode
            << ", priority = " << priority_
            << ", regions = " << regionIndices_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellSizeFunction> Foam::cellSizeFunction::New
(
    const dictionary& cellSizeFunctionDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize,
    const labelList regionIndices
)
{
    const word functionName
    (
        cellSizeFunctionDict.get<word>("cellSizeFunction")
    );

    Info<< indent << "Selecting cellSizeFunction "
        << functionName << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(functionName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown cellSizeFunction type "
            << functionName << nl << nl
            << "Valid cellSizeFunction types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<cellSizeFunction>
    (
        cstrIter()
        (
            cellSizeFunctionDict,
            surface,
            defaultCellSize,
            regionIndices
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeFunction::~cellSizeFunction()
{}


// ************************************************************************* //
