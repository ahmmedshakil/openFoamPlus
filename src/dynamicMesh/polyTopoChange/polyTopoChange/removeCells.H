/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::removeCells

Description
    Given list of cells to remove, insert all the topology changes.

    Works in two passes:
    - get faces that will become boundary faces
    - given these faces and the patches they have to go into make the
      changes.

SourceFiles
    removeCells.C

\*---------------------------------------------------------------------------*/

#ifndef removeCells_H
#define removeCells_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class bitSet;
class polyMesh;
class polyTopoChange;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class removeCells Declaration
\*---------------------------------------------------------------------------*/

class removeCells
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether or not to synchronize parallel case.
        const bool syncPar_;

public:

    //- Runtime type information
    ClassName("removeCells");


    // Constructors

        //- Construct from mesh. With parallel synchronization.
        explicit removeCells(const polyMesh& mesh);

        //- Construct from mesh, optionally with parallel synchronization.
        removeCells(const polyMesh& mesh, const bool syncPar);


    //- Destructor
    ~removeCells() = default;


    // Member Functions

        //- Get labels of faces exposed after cells removal.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the sides
        //    gets deleted)
        labelList getExposedFaces(const bitSet& removedCell) const;

        //- Get labels of faces exposed after cells removal.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the sides
        //    gets deleted)
        labelList getExposedFaces(const labelUList& cellsToRemove) const;

        //- Play commands into polyTopoChange to remove cells.
        //  patchIDs is for every element in facesToExpose (see above) the
        //  patch it has to go into. This cannot be a coupled patch!
        void setRefinement
        (
            const bitSet& removedCell,
            const labelUList& facesToExpose,
            const labelUList& patchIDs,
            polyTopoChange&
        ) const;

        //- Play commands into polyTopoChange to remove cells.
        //  patchIDs is for every element in facesToExpose (see above) the
        //  patch it has to go into. This cannot be a coupled patch!
        void setRefinement
        (
            const labelUList& cellsToRemove,
            const labelUList& facesToExpose,
            const labelUList& patchIDs,
            polyTopoChange&
        ) const;

        //- Force recalculation of locally stored data on topological change
        void updateMesh(const mapPolyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
