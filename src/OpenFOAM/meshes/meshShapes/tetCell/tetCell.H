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

Class
    Foam::tetCell

Description
    A tetrahedral cell primitive.

    It is important that the ordering of edges is the same for a tetrahedron
    class, a tetrahedron cell shape model and a tetCell

SourceFiles
    tetCell.C
    tetCellI.H

\*---------------------------------------------------------------------------*/

#ifndef tetCell_H
#define tetCell_H

#include "FixedList.H"
#include "triFace.H"
#include "edge.H"
#include "pointField.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellShape;

/*---------------------------------------------------------------------------*\
                           class tetCell Declaration
\*---------------------------------------------------------------------------*/

class tetCell
:
    public FixedList<label, 4>
{

public:

    // Constructors

        //- Construct null, with invalid point labels (-1)
        inline tetCell();

        //- Construct from four point labels
        inline tetCell
        (
            const label a,
            const label b,
            const label c,
            const label d
        );

        //- Construct from FixedList of four point labels
        inline tetCell(const FixedList<label, 4>& lst);

        //- Construct from an initializer list of four point labels
        inline explicit tetCell(std::initializer_list<label> lst);

        //- Construct from Istream
        inline tetCell(Istream& is);


    // Member Functions

        // Access

            //- Return i-th face
            inline triFace face(const label facei) const;

            //- Return first face adjacent to the given edge
            inline label edgeFace(const label edgei) const;

            //- Return face adjacent to the given face sharing the same edge
            inline label edgeAdjacentFace
            (
                const label edgei,
                const label facei
            ) const;

            //- Return i-th edge
            inline edge tetEdge(const label edgei) const;


        // Operations

            //- Return tet shape cell
            cellShape tetCellShape() const;

            //- Return the tetrahedron
            inline tetPointRef tet(const UList<point>& points) const;
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Data associated with the type are contiguous
template<>
inline bool contiguous<tetCell>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
