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

Class
    Foam::Vector2D

Description
    Templated 2D Vector derived from VectorSpace adding construction from
    2 components, element access using x() and y() member functions and
    the inner-product (dot-product).

SourceFiles
    Vector2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Vector2D_H
#define Vector2D_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Vector2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Vector2D
:
    public VectorSpace<Vector2D<Cmpt>, Cmpt, 2>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Vector2D<label> labelType;


    // Member constants

        //- Rank of Vector2D is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y };


    // Constructors

        //- Construct null
        inline Vector2D();

        //- Construct initialized to zero
        inline Vector2D(const Foam::zero);

        //- Copy construct from VectorSpace of the same rank
        inline Vector2D(const VectorSpace<Vector2D<Cmpt>, Cmpt, 2>& vs);

        //- Construct from two components
        inline Vector2D(const Cmpt& vx, const Cmpt& vy);

        //- Construct from Istream
        inline Vector2D(Istream& is);


    // Member Functions

        //- Access to the vector x component
        inline const Cmpt& x() const;

        //- Access to the vector y component
        inline const Cmpt& y() const;

        //- Access to the vector x component
        inline Cmpt& x();

        //- Access to the vector y component
        inline Cmpt& y();


        //- Normalise the vector by its magnitude
        inline Vector2D<Cmpt>& normalise();


        //- Perp dot product (dot product with perpendicular vector)
        inline scalar perp(const Vector2D<Cmpt>& b) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "Vector2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
