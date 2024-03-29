/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2008-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::tetherPotential> Foam::tetherPotential::New
(
    const word& name,
    const dictionary& propDict
)
{
    const word potentialType(propDict.get<word>("tetherPotential"));

    Info<< nl << "Selecting tether potential " << potentialType
        << " for " << name << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(potentialType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown tetherPotential type "
            << potentialType << nl << nl
            << "Valid tetherPotential types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<tetherPotential>(cstrIter()(name, propDict));
}


// ************************************************************************* //
