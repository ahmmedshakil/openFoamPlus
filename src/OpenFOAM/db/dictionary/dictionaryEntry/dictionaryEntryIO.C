/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Modified code Copyright (C) 2016 OpenCFD Ltd.
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

#include "dictionaryEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    entry(keyType(is)),
    dictionary(parentDict, is)
{
    is.fatalCheck(FUNCTION_NAME);
}


Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    Istream& is
)
:
    entry(key),
    dictionary(key, parentDict, is)
{
    is.fatalCheck(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dictionaryEntry::write(Ostream& os) const
{
    dictionary::writeEntry(keyword(), os);
}


// * * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const dictionaryEntry& e)
{
    e.write(os);
    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<dictionaryEntry>& ip
)
{
    const dictionaryEntry& e = ip.t_;

    os  << "    dictionaryEntry '" << e.keyword() << "'" << endl;

    return os;
}


// ************************************************************************* //
