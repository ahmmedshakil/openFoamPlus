/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
     Point addressing on the patch: pointEdges and pointFaces.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "SLList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcPointEdges() const
{
    DebugInFunction
        << "Calculating pointEdges" << endl;

    if (pointEdgesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "pointEdges already calculated"
            << abort(FatalError);
    }

    pointEdgesPtr_ = new labelListList(meshPoints().size());

    labelListList& pe = *pointEdgesPtr_;

    invertManyToMany(pe.size(), edges(), pe);

    DebugInfo
        << "    Finished." << endl;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcPointFaces() const
{
    DebugInFunction
        << "Calculating pointFaces" << endl;

    if (pointFacesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "pointFaces already calculated"
            << abort(FatalError);
    }

    const List<Face>& f = localFaces();

    // set up storage for pointFaces
    List<SLList<label>> pointFcs(meshPoints().size());

    forAll(f, facei)
    {
        const Face& curPoints = f[facei];

        for (const label pointi : curPoints)
        {
            pointFcs[pointi].append(facei);
        }
    }

    // sort out the list
    pointFacesPtr_ = new labelListList(pointFcs.size());

    labelListList& pf = *pointFacesPtr_;

    forAll(pointFcs, pointi)
    {
        pf[pointi].setSize(pointFcs[pointi].size());

        label i = 0;
        for (const label facei : pointFcs[pointi])
        {
            pf[pointi][i++] = facei;
        }
    }

    DebugInfo
        << "    Finished." << endl;
}


// ************************************************************************* //
