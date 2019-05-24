/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Typedef
    Foam::tensor2D

Description
    Tensor2D of scalars, i.e. Tensor2D<scalar>.

    Analytical computation of eigenvalues and eigenvectors of given tensor2D.

    Reference:
    \verbatim
        2-by-2 eigenvalue algorithm:
            Ford, W. (2014).
            Numerical linear algebra with applications: Using MATLAB.
            London: Elsevier/Academic Press.
            doi:10.1016/C2011-0-07533-6

            Blinn, J. (1996).
            Consider the lowly 2 x 2 matrix.
            IEEE Computer Graphics and Applications, 16(2), 82-88.
            doi:10.1109/38.486688

        2-by-2 eigenvector algorithm:
            Knill, O. (2004).
            Mathematics Math21b Fall 2004.
            bit.ly/2kjPVlX (Retrieved:06-09-19)
        
        Kahan summation algorithm for 2-by-2 matrix determinants:
            Jeannerod, C.-P., Louvet, N., & Muller, J.-M., (2013).
            Further analysis of Kahan's algorithm for the accurate computation
            of 2x2 determinants.
            Math. Comp. 82 (2013), 2245-2264.
            doi:10.1090/S0025-5718-2013-02679-8
    \endverbatim

SourceFiles
    tensor2D.C

\*---------------------------------------------------------------------------*/

#ifndef tensor2D_H
#define tensor2D_H

#include "Tensor2D.H"
#include "vector2D.H"
#include "complex.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Tensor2D<scalar> tensor2D;

//- Return eigenvalues of a given tensor2D
//  \param T the input tensor2D
//  \return Vector2D<complex> eigenvalues
Vector2D<complex> eigenValues
(
    const tensor2D& T
);


//- Return an eigenvector of a given tensor2D in a given direction1
//- corresponding to a given eigenvalue
//  \param T the input tensor2D
//  \param lambda the input eigenvalue
//  \param direction1 the input tensor2D component
//  \return Vector2D<complex> an eigenvector
Vector2D<complex> eigenVector
(
    const tensor2D& T,
    const complex& lambda,
    const Vector2D<complex>& direction
);


//- Return eigenvectors of a given tensor2D corresponding to given eigenvalues
//  \param T the input tensor2D
//  \param lambdas the input eigenvalues
//  \return Tensor2D<complex> eigenvectors (each row is an eigenvector)
Tensor2D<complex> eigenVectors
(
    const tensor2D& T,
    const Vector2D<complex>& lambdas
);


//- Return eigenvectors of a given tensor2D by
//- inherently computing the eigenvalues of the tensor2D
//  \param T the input tensor2D
//  \return Tensor2D<complex> eigenvectors (each row is an eigenvector)
Tensor2D<complex> eigenVectors
(
    const tensor2D& T
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sphericalTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
