/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "pointToFace.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToFace, word);
    addToRunTimeSelectionTable(topoSetSource, pointToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, pointToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, pointToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        pointToFace,
        word,
        point
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        pointToFace,
        istream,
        point
    );
}


Foam::topoSetSource::addToUsageTable Foam::pointToFace::usage_
(
    pointToFace::typeName,
    "\n    Usage: pointToFace <pointSet> any|all|edge\n\n"
    "    Select faces with\n"
    "    -any point in the pointSet\n"
    "    -all points in the pointSet\n\n"
    "    -two consecutive points (an edge) in the pointSet\n\n"
);

const Foam::Enum
<
    Foam::pointToFace::pointAction
>
Foam::pointToFace::pointActionNames_
({
    { pointAction::ANY, "any" },
    { pointAction::ALL, "all" },
    { pointAction::EDGE, "edge" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToFace::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName);

    const labelHashSet& pointLabels = loadedSet;

    if (option_ == ANY)
    {
        // Add faces with any point in loadedSet
        for (const label pointi : pointLabels)
        {
            const labelList& pFaces = mesh_.pointFaces()[pointi];

            for (const label facei : pFaces)
            {
                addOrDelete(set, facei, add);
            }
        }
    }
    else if (option_ == ALL)
    {
        // Add all faces whose points are all in set.

        // Count number of points using face.
        Map<label> numPoints(pointLabels.size());

        for (const label pointi : pointLabels)
        {
            const labelList& pFaces = mesh_.pointFaces()[pointi];

            for (const label facei : pFaces)
            {
                ++(numPoints(facei, 0));
            }
        }


        // Include faces that are referenced as many times as there are points
        // in face -> all points of face
        forAllConstIters(numPoints, iter)
        {
            const label facei = iter.key();
            const label count = iter.val();

            if (count == mesh_.faces()[facei].size())
            {
                addOrDelete(set, facei, add);
            }
        }
    }
    else if (option_ == EDGE)
    {
        const faceList& faces = mesh_.faces();

        forAll(faces, facei)
        {
            const face& f = faces[facei];

            forAll(f, fp)
            {
                if
                (
                    pointLabels.found(f[fp])
                 && pointLabels.found(f.nextLabel(fp))
                )
                {
                    addOrDelete(set, facei, add);
                    break;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetFaceSource(mesh),
    names_(one(), setName),
    option_(option)
{}


Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh),
    names_(),
    option_(pointActionNames_.get("option", dict))
{
    // Look for 'sets' or 'set'
    if (!dict.readIfPresent("sets", names_))
    {
        names_.resize(1);
        dict.readEntry("set", names_.first());
    }
}


Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    names_(one(), word(checkIs(is))),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces according to pointSet "
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, true, setName);
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing faces according to pointSet "
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, false, setName);
        }
    }
}


// ************************************************************************* //
