/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "IOstream.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::fileName Foam::IOstream::staticName_("IOstream");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::IOstream::name() const
{
    return staticName_;
}


Foam::fileName& Foam::IOstream::name()
{
    return staticName_;
}


bool Foam::IOstream::check(const char* operation) const
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "error in IOstream " << name() << " for operation " << operation
            << exit(FatalIOError);
    }

    return !bad();
}


void Foam::IOstream::fatalCheck(const char* operation) const
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "error in IOstream " << name() << " for operation " << operation
            << exit(FatalIOError);
    }
}


void Foam::IOstream::print(Ostream& os) const
{
    os  << "IOstream: " << "Version "  << version() << ", format "
        << format() << ", line " << lineNumber();

    if (opened())
    {
        os  << ", OPENED";
    }

    if (closed())
    {
        os  << ", CLOSED";
    }

    if (good())
    {
        os  << ", GOOD";
    }

    if (eof())
    {
        os  << ", EOF";
    }

    if (fail())
    {
        os  << ", FAIL";
    }

    if (bad())
    {
        os  << ", BAD";
    }

    os  << endl;
}


void Foam::IOstream::print(Ostream& os, const int streamState) const
{
    if (streamState == ios_base::goodbit)
    {
        os  << "ios_base::goodbit set : the last operation on stream succeeded"
            << endl;
    }
    else if (streamState & ios_base::badbit)
    {
        os  << "ios_base::badbit set : characters possibly lost"
            << endl;
    }
    else if (streamState & ios_base::failbit)
    {
        os  << "ios_base::failbit set : some type of formatting error"
            << endl;
    }
    else if (streamState & ios_base::eofbit)
    {
        os  << "ios_base::eofbit set : at end of stream"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<IOstream>& ip)
{
    ip.t_.print(os);
    return os;
}


// ************************************************************************* //
