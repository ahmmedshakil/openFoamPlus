/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::blockEdges::BSplineEdge

Description
    A blockEdge interface for B-splines.

SourceFiles
    BSplineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_BSplineEdge_H
#define blockEdges_BSplineEdge_H

#include "blockEdge.H"
#include "BSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                      Class BSplineEdge Declaration
\*---------------------------------------------------------------------------*/

class BSplineEdge
:
    public blockEdge,
    public BSpline
{
    // Private Member Functions

        //- No copy construct
        BSplineEdge(const BSplineEdge&) = delete;

        //- No copy assignment
        void operator=(const BSplineEdge&) = delete;


public:

    //- Runtime type information
    TypeName("BSpline");


    // Constructors

        //- Construct from components
        BSplineEdge
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& internalPoints
        );

        //- Construct from Istream, setting pointsList
        BSplineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField&,
            Istream&
        );


    //- Destructor
    virtual ~BSplineEdge() = default;


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- Return the length of the spline curve (not implemented)
        virtual scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
