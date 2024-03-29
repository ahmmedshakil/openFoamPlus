/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::isoSurfaceCell

Description
    A surface formed by the iso value.
    After "Polygonising A Scalar Field Using Tetrahedrons", Paul Bourke
    (http://paulbourke.net/geometry/polygonise) and
    "Regularised Marching Tetrahedra: improved iso-surface extraction",
    G.M. Treece, R.W. Prager and A.H. Gee.

    See isoSurface. This is a variant which does tetrahedrisation from
    triangulation of face to cell centre instead of edge of face to two
    neighbouring cell centres. This gives much lower quality triangles
    but they are local to a cell.

SourceFiles
    isoSurfaceCell.C
    isoSurfaceCellTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceCell_H
#define isoSurfaceCell_H

#include "labelPair.H"
#include "pointIndexHit.H"
#include "bitSet.H"
#include "boundBox.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class polyMesh;
class triSurface;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceCell Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceCell
:
    public meshedSurface
{
    // Private typedefs for convenience
    typedef meshedSurface MeshStorage;


    // Private data

        enum segmentCutType
        {
            NEARFIRST,      //!< Intersection close to e.first()
            NEARSECOND,     //!< Intersection close to e.second()
            NOTNEAR         //!< No intersection
        };

        enum cellCutType
        {
            NOTCUT,     //!< Cell not cut
            SPHERE,     //!< All edges to cell centre cut
            CUT         //!< Normal cut
        };


        //- Reference to mesh
        const polyMesh& mesh_;

        const scalarField& cVals_;

        const scalarField& pVals_;

        //- isoSurfaceCell value
        const scalar iso_;

        //- Optional bounds
        const boundBox bounds_;

        //- Optional cells to ignore.
        const bitSet& ignoreCells_;

        //- When to merge points
        const scalar mergeDistance_;

        //- Whether cell might be cut
        List<cellCutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- For every triangle the original cell in mesh
        labelList meshCells_;

        //- For every unmerged triangle point the point in the triSurface
        labelList triPointMergeMap_;

        //- triSurface points that have weighted interpolation
        DynamicList<label> interpolatedPoints_;

        //- corresponding original, unmerged points
        DynamicList<FixedList<label, 3>> interpolatedOldPoints_;

        //- corresponding weights
        DynamicList<FixedList<scalar, 3>> interpolationWeights_;


    // Private Member Functions

        //- Get location of iso value as fraction inbetween s0,s1
        scalar isoFraction
        (
            const scalar s0,
            const scalar s1
        ) const;

        //- Does any edge of triangle cross iso value?
        bool isTriCut
        (
            const triFace& tri,
            const scalarField& pointValues
        ) const;

        //- Determine whether cell is cut
        cellCutType calcCutType
        (
            const bitSet&,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const label
        ) const;

        void calcCutTypes
        (
            const bitSet&,
            const scalarField& cellValues,
            const scalarField& pointValues
        );

        //- Return the two common points between two triangles
        static labelPair findCommonPoints
        (
            const labelledTri&,
            const labelledTri&
        );

        //- Calculate centre of surface.
        static point calcCentre(const triSurface&);

        //- Replace surface (localPoints, localTris) with single point.
        //  Returns point.  Destroys arguments.
        pointIndexHit collapseSurface
        (
            const label celli,
            pointField& localPoints,
            DynamicList<labelledTri, 64>& localTris
        ) const;

        //- Determine per cc whether all near cuts can be snapped to single
        //  point.
        void calcSnappedCc
        (
            const bitSet& isTet,
            const scalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& triPoints,
            labelList& snappedCc
        ) const;

        //- Generate triangles for face connected to pointi
        void genPointTris
        (
            const scalarField& cellValues,
            const scalarField& pointValues,
            const label pointi,
            const label facei,
            const label celli,
            DynamicList<point, 64>& localTriPoints
        ) const;

        //- Generate triangles for tet connected to pointi
        void genPointTris
        (
            const scalarField& pointValues,
            const label pointi,
            const label facei,
            const label celli,
            DynamicList<point, 64>& localTriPoints
        ) const;

        //- Determine per point whether all near cuts can be snapped to single
        //  point.
        void calcSnappedPoint
        (
            const bitSet& isTet,
            const scalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& triPoints,
            labelList& snappedPoint
        ) const;

        //- Generate single point by interpolation or snapping
        template<class Type>
        Type generatePoint
        (
            const DynamicList<Type>& snappedPoints,
            const scalar s0,
            const Type& p0,
            const label p0Index,
            const scalar s1,
            const Type& p1,
            const label p1Index
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const DynamicList<Type>& snapped,
            const scalar s0,
            const Type& p0,
            const label p0Index,
            const scalar s1,
            const Type& p1,
            const label p1Index,
            const scalar s2,
            const Type& p2,
            const label p2Index,
            const scalar s3,
            const Type& p3,
            const label p3Index,
            DynamicList<Type>& points
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const scalarField& cVals,
            const scalarField& pVals,

            const Field<Type>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        triSurface stitchTriPoints
        (
            const bool checkDuplicates,
            const List<point>& triPoints,
            labelList& triPointReverseMap,  // unmerged to merged point
            labelList& triMap               // merged to unmerged triangle
        ) const;

        //- Determine edge-face addressing
        void calcAddressing
        (
            const triSurface& surf,
            List<FixedList<label, 3>>& faceEdges,
            labelList& edgeFace0,
            labelList& edgeFace1,
            Map<labelList>& edgeFacesRest
        ) const;

        //- Is triangle (given by 3 edges) not fully connected?
        static bool danglingTriangle
        (
            const FixedList<label, 3>& fEdges,
            const labelList& edgeFace1
        );

        //- Mark all non-fully connected triangles
        static label markDanglingTriangles
        (
            const List<FixedList<label, 3>>& faceEdges,
            const labelList& edgeFace0,
            const labelList& edgeFace1,
            const Map<labelList>& edgeFacesRest,
            boolList& keepTriangles
        );

        static triSurface subsetMesh
        (
            const triSurface&,
            const labelList& newToOldFaces,
            labelList& oldToNewPoints,
            labelList& newToOldPoints
        );


public:

    //- Runtime type information
    TypeName("isoSurfaceCell");


    // Constructors

        //- Construct from cell and point values
        //
        //  \param bounds optional bounding box for trimming
        //  \param mergeTol fraction of mesh bounding box for merging points
        //  \param ignoreCells cells to ignore in the cellValues
        isoSurfaceCell
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const bool regularise,
            const boundBox& bounds = boundBox::invertedBox,
            const scalar mergeTol = 1e-6,
            const bitSet& ignoreCells = bitSet()
        );


    // Member Functions

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            return meshCells_;
        }


        //- Interpolates cCoords, pCoords.
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceCellTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
