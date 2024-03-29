/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::hierarchGeomDecomp

Description
    Does hierarchical decomposition of points, selectable as \c hierarchical.

    Works by first sorting the points in x direction into equal sized bins,
    then in y direction and finally in z direction.

    Uses single array to hold decomposition which is indexed as if it is a
    3 dimensional array:

        finalDecomp[i,j,k] is indexed as

        i*n[0]*n[1] + j*n[1] + k

    E.g. if we're sorting 'xyz': the first sort (over the x-component)
    determines in which x-domain the point goes. Then for each of the x-domains
    the points are sorted in y direction and each individual x-domain gets
    split into three y-domains. And similar for the z-direction.

    Since the domains are of equal size the maximum difference in size is
    n[0]*n[1] (or n[1]*n[2]?) (small anyway)

    Method coefficients:
    \table
        Property  | Description                             | Required | Default
        n         | (nx ny nz)                              | yes |
        delta     | delta for rotation matrix               | no  | 0.001
        order     | order of operation                      | no  | xyz
    \endtable

SourceFiles
    hierarchGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef hierarchGeomDecomp_H
#define hierarchGeomDecomp_H

#include "geomDecomp.H"
#include "FixedList.H"
#include "direction.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hierarchGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class hierarchGeomDecomp
:
    public geomDecomp
{
    // Private Data

        //- Decomposition order in terms of components.
        FixedList<direction, 3> order_;


    // Private Member Functions

        //- Convert ordering string ("xyz") into list of components.
        //  Checks for bad entries, but no check for duplicate entries.
        void setOrder();

        //- Find index of value in list between
        //- first (inclusive) and last (exclusive)
        static label findLower
        (
            const UList<scalar>& list,
            const scalar val,
            const label first,
            const label last
        );

        //- Evaluates the weighted sizes for each sorted point.
        static void calculateSortedWeightedSizes
        (
            const labelList& current,
            const labelList& indices,
            const scalarField& weights,
            const label globalCurrentSize,

            scalarField& sortedWeightedSizes
        );

        //- Find midValue (at local index mid) such that the number of
        //  elements between mid and leftIndex are (globally summed) the
        //  wantedSize. Binary search.
        //
        //  \Return False if the binary search completed
        static bool findBinary
        (
            const label sizeTol,        // size difference considered acceptable
            const List<scalar>&,
            const label leftIndex,      // index of previous value
            const scalar leftValue,     // value at leftIndex
            const scalar maxValue,      // global max of values
            const scalar wantedSize,    // wanted size
            label& mid,                 // index where size of bin is wantedSize
            scalar& midValue            // value at mid
        );

        //- Find midValue (at local index mid) such that the number of
        //  elements between mid and leftIndex are (globally summed) the
        //  wantedSize. Binary search.
        static bool findBinary
        (
            const label sizeTol,        // size difference considered acceptable
            const List<scalar>& sortedWeightedSizes,
            const List<scalar>&,
            const label leftIndex,      // index of previous value
            const scalar leftValue,     // value at leftIndex
            const scalar maxValue,      // global max of values
            const scalar wantedSize,    // wanted size
            label& mid,                 // index where size of bin is wantedSize
            scalar& midValue            // value at mid
        );

        //- Recursively sort in x,y,z (or rather acc. to decompOrder_)
        //  \return the number of warnings from findBinary
        label sortComponent
        (
            const label sizeTol,
            const pointField&,
            const labelList& slice,         // slice of points to decompose
            const direction componentIndex, // index in decompOrder_
            const label prevMult,           // multiplication factor
            labelList& finalDecomp          // overall decomposition
        ) const;

        //- Recursively sort in x,y,z (or rather acc. to decompOrder_)
        //- Using weighted points.
        //  \return the number of warnings from findBinary
        label sortComponent
        (
            const label sizeTol,
            const scalarField& weights,
            const pointField&,
            const labelList& slice,         // slice of points to decompose
            const direction componentIndex, // index in decompOrder_
            const label prevMult,           // multiplication factor
            labelList& finalDecomp          // overall decomposition
        ) const;


        //- No copy construct
        hierarchGeomDecomp(const hierarchGeomDecomp&) = delete;

        //- No copy assignment
        void operator=(const hierarchGeomDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("hierarchical");


    // Constructors

        //- Construct given the decomposition dictionary
        hierarchGeomDecomp(const dictionary& decompDict);

        //- Construct for decomposition dictionary and region name
        hierarchGeomDecomp
        (
            const dictionary& decompDict,
            const word& regionName
        );


    //- Destructor
    virtual ~hierarchGeomDecomp() = default;


    // Member Functions

        //- Hierarchical is aware of processor boundaries
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField&,
            const scalarField& weights
        ) const;

        //- Decompose with uniform weights.
        //  Code for weighted decomposition is a bit complex,
        //  so kept separate for now.
        virtual labelList decompose(const pointField&) const;


        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed).
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights
        ) const
        {
            checkDecompositionDirections(mesh.geometricD());
            return decompose(cc, cWeights);
        }

        //- Decompose with uniform weights.
        //  Code for weighted decomposition is a bit complex,
        //  so kept separate for now.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const
        {
            checkDecompositionDirections(mesh.geometricD());
            return decompose(cc);
        }

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        ) const
        {
            return decompose(cc, cWeights);
        }

        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc
        ) const
        {
            return decompose(cc);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
