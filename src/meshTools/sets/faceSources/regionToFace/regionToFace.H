/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::regionToFace

Description
    A topoSetFaceSource to select faces belonging to a topological connected
    region (that contains given point)

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        set         | The face set name to use          | yes   |
        nearPoint   | The point on/near to the region   | yes   |
    \endtable

SourceFiles
    regionToFace.C

\*---------------------------------------------------------------------------*/

#ifndef regionToFace_H
#define regionToFace_H

#include "topoSetFaceSource.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class regionToFace Declaration
\*---------------------------------------------------------------------------*/

class regionToFace
:
    public topoSetFaceSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

        //- Coordinate that is nearest/on connected region
        point nearPoint_;


    // Private Member Functions

        //- Walk edge-face-edge
        void markZone
        (
            const indirectPrimitivePatch& patch,
            const label proci,
            const label facei,
            const label zoneI,
            labelList& faceZone
        ) const;

        void combine(topoSet& set, const bool add) const;

public:

    //- Runtime type information
    TypeName("regionToFace");

    // Constructors

        //- Construct from components
        regionToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const point& nearPoint
        );

        //- Construct from dictionary
        regionToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        regionToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~regionToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
