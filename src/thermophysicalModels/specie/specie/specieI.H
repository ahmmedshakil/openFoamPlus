/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline specie::specie
(
    const word& name,
    const scalar Y,
    const scalar molWeight
)
:
    name_(name),
    Y_(Y),
    molWeight_(molWeight)
{}


inline specie::specie
(
    const scalar Y,
    const scalar molWeight
)
:
    Y_(Y),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline specie::specie(const word& name, const specie& st)
:
    name_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& specie::name() const
{
    return name_;
}


inline scalar specie::W() const
{
    return molWeight_;
}


inline scalar specie::Y() const
{
    return Y_;
}


inline scalar specie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void specie::operator=(const specie& st)
{
    //name_ = st.name_;
    Y_ = st.Y_;
    molWeight_ = st.molWeight_;
}


inline void specie::operator+=(const specie& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > SMALL)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumY;
}


inline void specie::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline specie operator+(const specie& st1, const specie& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, SMALL);

    if (mag(sumY) > SMALL)
    {
        return specie
        (
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_)
        );
    }
    else
    {
        return st1;
    }
}


inline specie operator*(const scalar s, const specie& st)
{
    return specie
    (
        s*st.Y_,
        st.molWeight_
    );
}


inline specie operator==(const specie& st1, const specie& st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < SMALL)
    {
        diffY = SMALL;
    }

    const scalar diffRW = st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;

    #ifdef __clang__
    // Using intermediate volatile bool to prevent compiler optimising out the
    // if block (above) - CLANG 3.7.1
    volatile const bool valid = (mag(diffRW) > SMALL);
    const scalar molWeight = valid ? diffY/diffRW : GREAT;
    #else
    scalar molWeight = GREAT;
    if (mag(diffRW) > SMALL)
    {
         molWeight = diffY/diffRW;
    }
    #endif

    return specie(diffY, molWeight);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
