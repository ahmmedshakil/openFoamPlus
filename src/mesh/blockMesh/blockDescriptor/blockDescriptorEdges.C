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

\*---------------------------------------------------------------------------*/

#include "blockDescriptor.H"
#include "lineEdge.H"
#include "lineDivide.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::blockDescriptor::edgePointsWeights
(
    pointField (&edgePoints)[12],
    scalarList (&edgeWeights)[12],
    const label edgei,
    const label start,
    const label end,
    const label nDiv
) const
{
    // Set reference to the list of labels defining the block
    const labelList& blockLabels = blockShape_;

    // Get list of points for this block
    const pointField blockPoints = blockShape_.points(vertices_);

    // Set the edge points/weights
    // The edge is a straight-line if it is not in the list of blockEdges

    for (const blockEdge& cedge : blockEdges_)
    {
        const int cmp = cedge.compare(blockLabels[start], blockLabels[end]);

        if (cmp)
        {
            if (cmp > 0)
            {
                // Curve has the same orientation

                // Divide the line
                const lineDivide divEdge(cedge, nDiv, expand_[edgei]);

                edgePoints[edgei]  = divEdge.points();
                edgeWeights[edgei] = divEdge.lambdaDivisions();
            }
            else
            {
                // Curve has the opposite orientation

                // Divide the line
                const lineDivide divEdge(cedge, nDiv, expand_[edgei].inv());

                const pointField& p = divEdge.points();
                const scalarList& d = divEdge.lambdaDivisions();

                edgePoints[edgei].setSize(p.size());
                edgeWeights[edgei].setSize(d.size());

                label pn = p.size() - 1;
                forAll(p, pi)
                {
                    edgePoints[edgei][pi]  = p[pn - pi];
                    edgeWeights[edgei][pi] = 1 - d[pn - pi];
                }
            }

            // Found curved-edge: done
            return 1;
        }
    }


    // Not curved-edge: divide the edge as a straight line
    lineDivide divEdge
    (
        blockEdges::lineEdge(blockPoints, start, end),
        nDiv,
        expand_[edgei]
    );

    edgePoints[edgei]  = divEdge.points();
    edgeWeights[edgei] = divEdge.lambdaDivisions();

    return 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::blockDescriptor::edgesPointsWeights
(
    pointField (&edgePoints)[12],
    scalarList (&edgeWeights)[12]
) const
{
    const label ni = sizes().x();
    const label nj = sizes().y();
    const label nk = sizes().z();

    label nCurvedEdges = 0;

    // X-direction
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 0,  0, 1, ni);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 1,  3, 2, ni);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 2,  7, 6, ni);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 3,  4, 5, ni);

    // Y-direction
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 4,  0, 3, nj);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 5,  1, 2, nj);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 6,  5, 6, nj);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 7,  4, 7, nj);

    // Z-direction
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 8,  0, 4, nk);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 9,  1, 5, nk);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 10, 2, 6, nk);
    nCurvedEdges += edgePointsWeights(edgePoints, edgeWeights, 11, 3, 7, nk);

    return nCurvedEdges;
}


// ************************************************************************* //
