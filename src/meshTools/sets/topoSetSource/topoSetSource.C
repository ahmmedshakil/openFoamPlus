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

\*---------------------------------------------------------------------------*/

#include "topoSetSource.H"
#include "polyMesh.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetSource, 0);
    defineRunTimeSelectionTable(topoSetSource, word);
    defineRunTimeSelectionTable(topoSetSource, istream);
}


Foam::HashTable<Foam::string>* Foam::topoSetSource::usageTablePtr_ = nullptr;


const Foam::Enum
<
    Foam::topoSetSource::setAction
>
Foam::topoSetSource::actionNames
({
    { setAction::ADD, "add" },
    { setAction::SUBTRACT, "subtract" },
    { setAction::SUBSET, "subset" },
    { setAction::INVERT, "invert" },
    { setAction::CLEAR, "clear" },
    { setAction::NEW, "new" },
    { setAction::REMOVE, "remove" },
    { setAction::LIST, "list" },
    { setAction::SUBTRACT, "delete" },   // Compat (1806)
});


const Foam::string Foam::topoSetSource::illegalSource_
(
    "Illegal topoSetSource name"
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::topoSetSource::check(labelList& list, const label maxLabel)
{
    const label len = list.size();

    label nGood = 0;

    for (label i=0; i < len; ++i)
    {
        const label val = list[i];

        if (val >= 0 && val < maxLabel)
        {
            if (nGood != i)
            {
                list[nGood] = val;
            }
            ++nGood;
        }
    }

    const label nReject = (len - nGood);

    if (nReject)
    {
        list.resize(nGood);

        // Report?
    }

    return !nReject;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    auto cstrIter = wordConstructorTablePtr_->cfind(topoSetSourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type "
            << topoSetSourceType << nl << nl
            << "Valid topoSetSource types :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, dict));
}


Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    auto cstrIter = istreamConstructorTablePtr_->cfind(topoSetSourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type "
            << topoSetSourceType << nl << nl
            << "Valid topoSetSource types :" << endl
            << istreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, is));
}


Foam::Istream& Foam::topoSetSource::checkIs(Istream& is)
{
    if (!is.good() || is.eof())
    {
        FatalErrorInFunction
            << exit(FatalError);
    }

    return is;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const label id,
    const bool add
) const
{
    if (add)
    {
        set.set(id);
    }
    else
    {
        set.unset(id);
    }
}


void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const labelUList& labels,
    const bool add
) const
{
    if (add)
    {
        set.set(labels);
    }
    else
    {
        set.unset(labels);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetSource::topoSetSource(const polyMesh& mesh)
:
    mesh_(mesh),
    verbose_(true)
{}


// ************************************************************************* //
