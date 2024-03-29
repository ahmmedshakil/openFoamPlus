/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::regionToCell

Description
    A topoSetCellSource to select cells belonging to a topological connected
    region (that contains given points)

Usage
    In dictionary input:
    \verbatim
        // optional name of cellSet delimiting search
        set         c0;

        //- Number of cell layers to erode mesh to detect holes in the mesh
        //  Set to 0 if not used.
        nErode 3;

        // points inside region to select
        insidePoints ((1 2 3));
    \endverbatim

    Dictionary parameters:
    \table
        Property    | Description                           | Required | Default
        insidePoints | Points inside regions                | yes   |
        nErode      | Cell layers to erode to detect holes  | no | 0
        set         | Restrict to named cellSet             | no    | ""
    \endtable

SourceFiles
    regionToCell.C

\*---------------------------------------------------------------------------*/

#ifndef regionToCell_H
#define regionToCell_H

#include "topoSetCellSource.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

/*---------------------------------------------------------------------------*\
                        Class regionToCell Declaration
\*---------------------------------------------------------------------------*/

class regionToCell
:
    public topoSetCellSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cellSet to keep to
        const word setName_;

        //- Coordinate(s) that is inside connected region
        const pointField insidePoints_;

        //- Number of layers to erode
        const label nErode_;


    // Private Member Functions

        //- Mark faces inbetween selected and unselected elements
        void markRegionFaces
        (
            const boolList& selectedCell,
            boolList& regionFace
        ) const;

        //- Determine for every disconnected region in the mesh whether
        //  it contains a locationInMesh
        boolList findRegions(const bool verbose, const regionSplit&) const;

        //- Unselect regions not containing a locationInMesh
        void unselectOutsideRegions(boolList& selectedCell) const;

        //- Unselect one layer of cells from selectedCell
        void shrinkRegions(boolList& selectedCell) const;

        //- Erode a given number of layers from selectedCell. Remove any
        //  region that gets disconnected that way.
        void erode(boolList& selectedCell) const;

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("regionToCell");

    // Constructors

        //- Construct from components
        regionToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const pointField& insidePoints,
            const label nErode
        );

        //- Construct from dictionary
        regionToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        regionToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~regionToCell() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
