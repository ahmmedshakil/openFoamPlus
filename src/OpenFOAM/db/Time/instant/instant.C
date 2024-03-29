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

#include "instant.H"
#include "Time.H"
#include <cstdlib>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::instant::typeName = "instant";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instant::instant(scalar timeValue)
:
    Instant<word>(timeValue, Time::timeName(timeValue))
{}


Foam::instant::instant(const word& timeName)
:
    Instant<word>(std::atof(timeName.c_str()), timeName)
{}


Foam::instant::instant(word&& timeName)
:
    Instant<word>(std::atof(timeName.c_str()), std::move(timeName))
{}


// ************************************************************************* //
