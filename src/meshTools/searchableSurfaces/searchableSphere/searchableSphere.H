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
    Foam::searchableSphere

Description
    Searching on sphere

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        type        | sphere / searchableSphere             | selector |
        origin      | The origin (centre) of the sphere     | yes   |
        radius      | The (outside) radius of sphere        | yes   |
        centre      | Alternative for 'origin'              | no    |
    \endtable

SourceFiles
    searchableSphere.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSphere_H
#define searchableSphere_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class searchableSphere Declaration
\*---------------------------------------------------------------------------*/

class searchableSphere
:
    public searchableSurface
{
private:

    // Private Member Data

        //- Centre point of the sphere
        const point origin_;

        //- The outer radius of the sphere
        const scalar radius_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Find nearest point on sphere.
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        //- Find intersection with sphere
        void findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;


        //- No copy construct
        searchableSphere(const searchableSphere&) = delete;

        //- No copy assignment
        void operator=(const searchableSphere&) = delete;


public:

    //- Runtime type information
    TypeName("searchableSphere");


    // Constructors

        //- Construct from components
        searchableSphere
        (
            const IOobject& io,
            const point& centre,
            const scalar radius
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableSphere
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSphere() = default;


    // Member Functions

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether supports volume type (below)
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- What is type of points outside bounds
        virtual volumeType outsideVolumeType() const
        {
            return volumeType::OUTSIDE;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const
        {
            return tmp<pointField>::New(1, origin_);
        }

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const
        {
            return coordinates();
        }

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

            //- Determine type (inside/outside/mixed) for point.
            //  Unknown if cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
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
