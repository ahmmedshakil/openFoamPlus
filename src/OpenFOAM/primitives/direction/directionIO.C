/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "direction.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::direction Foam::readDirection(Istream& is)
{
    direction val(0);
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, direction& d)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get direction"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        d = direction(t.labelToken());
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (direction), found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const direction d)
{
    os.write(label(d));
    os.check(FUNCTION_NAME);
    return os;
}


std::ostream& Foam::operator<<(std::ostream& os, const direction d)
{
    os << int(d);
    return os;
}


// ************************************************************************* //
