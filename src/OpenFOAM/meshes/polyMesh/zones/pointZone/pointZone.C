/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "pointZone.H"
#include "addToRunTimeSelectionTable.H"
#include "pointZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointZone, 0);
    defineRunTimeSelectionTable(pointZone, dictionary);
    addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}

const char* const Foam::pointZone::labelsName = "pointLabels";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointZone::pointZone
(
    const word& name,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const word& name,
    const labelUList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const word& name,
    labelList&& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, std::move(addr), index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, dict, this->labelsName, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const pointZone& origZone,
    const labelUList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(origZone, addr, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const pointZone& origZone,
    labelList&& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(origZone, std::move(addr), index),
    zoneMesh_(zm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointZoneMesh& Foam::pointZone::zoneMesh() const
{
    return zoneMesh_;
}


Foam::label Foam::pointZone::whichPoint(const label globalPointID) const
{
    return zone::localID(globalPointID);
}


bool Foam::pointZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh_.mesh().points().size(), report);
}


bool Foam::pointZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zoneMesh().mesh();

    labelList maxZone(mesh.nPoints(), -1);
    labelList minZone(mesh.nPoints(), labelMax);

    const labelList& addr = *this;

    for (const label pointi : addr)
    {
        maxZone[pointi] = index();
        minZone[pointi] = index();
    }
    syncTools::syncPointList(mesh, maxZone, maxEqOp<label>(), label(-1));
    syncTools::syncPointList(mesh, minZone, minEqOp<label>(), labelMax);

    bool error = false;

    forAll(maxZone, pointi)
    {
        // Check point in same (or no) zone on all processors
        if
        (
            (
                maxZone[pointi] != -1
             || minZone[pointi] != labelMax
            )
         && (maxZone[pointi] != minZone[pointi])
        )
        {
            if (report && !error)
            {
                Info<< " ***Problem with pointZone " << index()
                    << " named " << name()
                    << ". Point " << pointi
                    << " at " << mesh.points()[pointi]
                    << " is in zone "
                    << (minZone[pointi] == labelMax ? -1 : minZone[pointi])
                    << " on some processors and in zone "
                    << maxZone[pointi]
                    << " on some other processors." << nl
                    << "(suppressing further warnings)"
                    << endl;
            }
            error = true;
        }
    }

    return error;
}


void Foam::pointZone::writeDict(Ostream& os) const
{
    os  << nl << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry(this->labelsName, os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pointZone::operator=(const pointZone& zn)
{
    clearAddressing();
    labelList::operator=(zn);
}


void Foam::pointZone::operator=(const labelUList& addr)
{
    clearAddressing();
    labelList::operator=(addr);
}


void Foam::pointZone::operator=(labelList&& addr)
{
    clearAddressing();
    labelList::transfer(addr);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointZone& zn)
{
    zn.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
