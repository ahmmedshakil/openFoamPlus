/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{
    defineTypeName(cartesian);
    addToRunTimeSelectionTable(coordinateSystem, cartesian, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSystem::cartesian::cartesian()
:
    coordinateSystem()
{}


Foam::coordSystem::cartesian::cartesian(const coordinateSystem& csys)
:
    coordinateSystem(csys)
{}


Foam::coordSystem::cartesian::cartesian(coordinateSystem&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::cartesian::cartesian(autoPtr<coordinateSystem>&& csys)
:
    coordinateSystem(std::move(csys))
{}


Foam::coordSystem::cartesian::cartesian
(
    const point& origin,
    const coordinateRotation& crot
)
:
    coordinateSystem(origin, crot)
{}


Foam::coordSystem::cartesian::cartesian
(
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(origin, axis, dirn)
{}


Foam::coordSystem::cartesian::cartesian
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::coordSystem::cartesian::cartesian
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


Foam::coordSystem::cartesian::cartesian(const dictionary& dict)
:
    coordinateSystem(dict)
{}


Foam::coordSystem::cartesian::cartesian
(
    const dictionary& dict,
    const word& dictName
)
:
    coordinateSystem(dict, dictName)
{}


// ************************************************************************* //
