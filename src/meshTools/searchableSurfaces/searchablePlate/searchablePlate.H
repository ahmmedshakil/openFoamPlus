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

Class
    Foam::searchablePlate

Description
    Searching on finite plate. Plate has to be aligned with coordinate
    axes.
    Plate defined as origin and span. One of the components of span has
    to be 0 which defines the normal direction. E.g.

    \verbatim
    span    = (Sx Sy 0)     // plate in x-y plane
    origin  = (Ox Oy Oz)
    \endverbatim

    now plane is from (Ox Oy Oz) to (Ox+Sx Oy+Sy Oz)

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        type        | plate / searchablePlate           | selector |
        origin      | centre of the plate               | yes   |
        span        | The plate dimensions              | yes   |
    \endtable

SourceFiles
    searchablePlate.C

\*---------------------------------------------------------------------------*/

#ifndef searchablePlate_H
#define searchablePlate_H

#include "searchableSurface.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class searchablePlate Declaration
\*---------------------------------------------------------------------------*/

class searchablePlate
:
    public searchableSurface
{
private:

    // Private Member Data

        const point origin_;

        const vector span_;

        //- Coordinate direction which is normal
        const direction normalDir_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Calculate normal direction from span
        static direction calcNormal(const point&);

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        pointIndexHit findLine
        (
            const point& start,
            const point& end
        ) const;

        //- No copy construct
        searchablePlate(const searchablePlate&) = delete;

        //- No copy assignment
        void operator=(const searchablePlate&) = delete;


public:

    //- Runtime type information
    TypeName("searchablePlate");


    // Constructors

        //- Construct from components
        searchablePlate
        (
            const IOobject& io,
            const point& origin,
            const vector& span
        );

        //- Construct from dictionary (used by searchableSurface)
        searchablePlate
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchablePlate() = default;


    // Member Functions

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return false;
        }

        //- What is type of points outside bounds
        virtual volumeType outsideVolumeType() const
        {
            return volumeType::UNKNOWN;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const;


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
