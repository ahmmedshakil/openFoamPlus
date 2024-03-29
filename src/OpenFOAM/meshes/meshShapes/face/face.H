/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::face

Description
    A face is a list of labels corresponding to mesh vertices.

See also
    Foam::triFace

SourceFiles
    faceI.H
    face.C
    faceIntersection.C
    faceContactSphere.C
    faceAreaInContact.C
    faceTemplates.C

\*---------------------------------------------------------------------------*/
#ifndef face_H
#define face_H

#include "pointField.H"
#include "labelList.H"
#include "edgeList.H"
#include "vectorField.H"
#include "faceListFwd.H"
#include "intersection.H"
#include "pointHit.H"
#include "FixedList.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class face;
class triFace;

template<class T, int SizeMin> class DynamicList;

inline Istream& operator>>(Istream& is, face& f);

/*---------------------------------------------------------------------------*\
                           Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public labelList
{
    // Private Member Functions

        //- Edge to the right of face vertex i
        inline label right(const label i) const;

        //- Edge to the left of face vertex i
        inline label left(const label i) const;

        //- Construct list of edge vectors for face
        tmp<vectorField> calcEdges
        (
            const UList<point>& points
        ) const;

        //- Cos between neighbouring edges
        scalar edgeCos
        (
            const vectorField& edges,
            const label index
        ) const;

        //- Find index of largest internal angle on face
        label mostConcaveAngle
        (
            const UList<point>& points,
            const vectorField& edges,
            scalar& edgeCos
        ) const;

        //- Enumeration listing the modes for split()
        enum splitMode
        {
            COUNTTRIANGLE,  //!< count if split into triangles
            COUNTQUAD,      //!< count if split into triangles and quads
            SPLITTRIANGLE,  //!< split into triangles
            SPLITQUAD       //!< split into  triangles and quads
        };

        //- Split face into triangles or triangles and quads.
        //  Stores results quadFaces[quadI], triFaces[triI]
        //  \return number of new faces created
        label split
        (
            const splitMode mode,
            const UList<point>& points,
            label& triI,
            label& quadI,
            faceList& triFaces,
            faceList& quadFaces
        ) const;


public:

    //- Return types for classify
    enum proxType
    {
        NONE,
        POINT,  // Close to point
        EDGE    // Close to edge
    };

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        inline face();

        //- Construct given size, with invalid point labels (-1)
        inline explicit face(const label sz);

        //- Copy construct from list of labels
        inline explicit face(const labelUList& list);

        //- Copy construct from list of labels
        template<unsigned N>
        inline explicit face(const FixedList<label, N>& list);

        //- Copy construct from an initializer list of labels
        inline explicit face(std::initializer_list<label> list);

        //- Move construct from list of labels
        inline explicit face(labelList&& list);

        //- Copy construct from triFace
        face(const triFace& f);

        //- Construct from Istream
        inline face(Istream& is);


    // Member Functions

        //- Collapse face by removing duplicate point labels
        //  return the collapsed size
        label collapse();

        //- Flip the face in-place.
        //  The starting points of the original and flipped face are identical.
        void flip();

        //- Return the points corresponding to this face
        inline pointField points(const UList<point>& points) const;

        //- Centre point of face
        point centre(const UList<point>& points) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average
        (
            const UList<point>& meshPoints,
            const Field<Type>& fld
        ) const;

        //- The area normal - with magnitude equal to area of face
        vector areaNormal(const UList<point>& p) const;

        //- The unit normal
        inline vector unitNormal(const UList<point>& p) const;

        //- Legacy name for areaNormal()
        //  \deprecated(2018-06) Deprecated for new use
        vector FOAM_DEPRECATED_FOR(2018-12, "areaNormal() or unitNormal()")
        normal(const UList<point>& p) const
        {
            return areaNormal(p); // Legacy definition
        }

        //- Magnitude of face area
        inline scalar mag(const UList<point>& p) const;

        //- Return face with reverse direction
        //  The starting points of the original and reverse face are identical.
        face reverseFace() const;

        // Navigation through face vertices

        //- Return true if the point label is found in face.
        inline bool found(const label pointLabel) const;

        //- Find local index on face for the point label,
        //  \return position in face (0,1,2,...) or -1 if not found.
        inline label which(const label pointLabel) const;

        //- Next vertex on face
        inline label nextLabel(const label i) const;

        //- Previous vertex on face
        inline label prevLabel(const label i) const;


        //- Return the volume swept out by the face when its points move
        scalar sweptVol
        (
            const UList<point>& oldPoints,
            const UList<point>& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //  point and density specification
        tensor inertia
        (
            const UList<point>& p,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return potential intersection with face with a ray starting
        //  at p, direction n (does not need to be normalized)
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. Face-centre is calculated from point average.
        //  For a hit, the distance is signed.  Positive number
        //  represents the point in front of triangle
        //  In case of miss the point is the nearest point on the face
        //  and the distance is the distance between the intersection point
        //  and the original point.
        //  The half-ray or full-ray intersection and the contact
        //  sphere adjustment of the projection vector is set by the
        //  intersection parameters
        pointHit ray
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. See triangle::intersection for details.
        pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const UList<point>& meshPoints,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        pointHit nearestPoint
        (
            const point& p,
            const UList<point>& meshPoints
        ) const;

        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        pointHit nearestPointClassify
        (
            const point& p,
            const UList<point>& meshPoints,
            label& nearType,
            label& nearLabel
        ) const;

        //- The sign for the side of the face plane the point is on,
        //- using three evenly distributed face points for the estimated normal.
        //  Uses the supplied tolerance for rounding around zero.
        //  \return
        //  -  0: on plane
        //  - +1: front-side
        //  - -1: back-side
        int sign
        (
            const point& p,
            const UList<point>& points,
            const scalar tol = SMALL
        ) const;

        //- Return contact sphere diameter
        scalar contactSphereDiameter
        (
            const point& p,
            const vector& n,
            const UList<point>& meshPoints
        ) const;

        //- Return area in contact, given the displacement in vertices
        scalar areaInContact
        (
            const UList<point>& meshPoints,
            const scalarField& v
        ) const;

        //- Return number of edges
        inline label nEdges() const;

        //- Return edges in face point ordering,
        //  i.e. edges()[0] is edge between [0] and [1]
        edgeList edges() const;

        //- Return n-th face edge
        inline edge faceEdge(const label n) const;

        //- The edge direction on the face
        //  \return
        //  -  0: edge not found on the face
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        int edgeDirection(const edge& e) const;

        //- Find the longest edge on a face.
        label longestEdge(const UList<point>& pts) const;

        // Face splitting utilities

            //- Number of triangles after splitting
            inline label nTriangles() const;

            //- Number of triangles after splitting
            label nTriangles(const UList<point>& unused) const;

            //- Split into triangles using existing points.
            //  Result in triFaces[triI..triI+nTri].
            //  Splits intelligently to maximize triangle quality.
            //  Returns number of faces created.
            label triangles
            (
                const UList<point>& points,
                label& triI,
                faceList& triFaces
            ) const;

            //- Split into triangles using existing points.
            //  Append to DynamicList.
            //  Returns number of faces created.
            template<int SizeMin>
            label triangles
            (
                const UList<point>& points,
                DynamicList<face, SizeMin>& triFaces
            ) const;

            //- Number of triangles and quads after splitting
            //  Returns the sum of both
            label nTrianglesQuads
            (
                const UList<point>& points,
                label& nTris,
                label& nQuads
            ) const;

            //- Split into triangles and quads.
            //  Results in triFaces (starting at triI) and quadFaces
            //  (starting at quadI).
            //  Returns number of new faces created.
            label trianglesQuads
            (
                const UList<point>& points,
                label& triI,
                label& quadI,
                faceList& triFaces,
                faceList& quadFaces
            ) const;

        //- Compare faces
        //  \return
        //  -  0: different
        //  - +1: identical
        //  - -1: same face, but different orientation
        static int compare(const face& a, const face& b);

        //- Return true if the faces have the same vertices
        static bool sameVertices(const face& a, const face& b);


    // Istream Operator

        friend Istream& operator>>(Istream& is, face& f);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Hash specialization for face
template<>
struct Hash<face>
{
    inline unsigned operator()(const face& obj, unsigned seed=0) const
    {
        return Hasher(obj.cdata(), obj.size()*sizeof(label), seed);
    }
};


//- Specialization to offset faces, used in ListListOps::combineOffset
template<>
struct offsetOp<face>
{
    inline face operator()
    (
        const face& x,
        const label offset
    ) const
    {
        face result(x.size());

        forAll(x, i)
        {
            result[i] = x[i] + offset;
        }
        return result;
    }
};


//- Deprecated(2017-04) find the longest edge on a face.
//- Face point labels index into pts.
//  \deprecated(2017-04) use class method instead
label longestEdge(const face& f, const UList<point>& pts);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const face& a, const face& b);
inline bool operator!=(const face& a, const face& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceI.H"

#ifdef NoRepository
    #include "faceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
