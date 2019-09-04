/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "vector2D.H"
#include "tensor2D.H"
#include "IFstream.H"
#include "complex.H"
#include <cmath>
#include "Random.H"
#include "RectangularMatrix.H"
#include "MatrixTools.H"

using namespace Foam;
#define RUNALL true
const bool verbose = true;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void horizontalLine()
{
    Info<< "+---------+---------+---------+---------+---------+" << nl;
}


// Create random tensor2D
tensor2D randomTensor2D
(
    Random& rnd
)
{
    tensor2D A(0.0);
    std::generate(A.begin(), A.end(), [&]{ return rnd.GaussNormal<scalar>(); });
    return A;
}


// Return true if two scalars are equal within a given tolerance, and v.v.
bool isEqual
(
    const scalar a,
    const scalar b,
    const bool verbose = false,
    const scalar relTol = 1e-5,
    const scalar absTol = 1e-8
)
{
    if ((absTol + relTol*mag(b)) < mag(a - b))
    {
        if (verbose)
        {
            Info<< "Elements are not close in terms of tolerances:"
                << nl << a << tab << b << nl;
        }
        return false;
    }

    if (verbose)
    {
        Info<< "All elems are the same within the tolerances" << nl;
    }

    return true;
}


// Check 2-by-2 eigen decomposition functions for tensor2D
void verification_eigendecomposition_tensor2D
(
    const tensor2D& T
)
{
    const Vector2D<complex> EVals(eigenValues(T));

    Info<< "# Eigenvalues of Tensor2D<scalar>: " << nl << EVals << nl;

    {
        Info<< "# Prod(Eigenvalues) = det(T):" << nl;
        // If the matrix entries are all real, then so is the determinant.
        const scalar determinant = det(T);
        complex EValsProd = pTraits<complex>::one;
        for (const auto& x : EVals)
        {
            EValsProd *= x;
        }
        isEqual(EValsProd.real(), determinant, verbose);
    }

    {
        Info<< "# Sum(Eigenvalues) = trace(T):" << nl;
        const scalar trace = T.xx() + T.yy();
        complex EValsSum = Zero;
        for (const auto& val : EVals)
        {
            EValsSum += val;
        }
        isEqual(EValsSum.real(), trace, verbose);
    }

    const Tensor2D<complex> EVecs(eigenVectors(T));

    Info<< "# Eigenvectors of Tensor2D<scalar>: " << nl << EVecs << nl;

    {
        Info<< "# T*eigenvector - eigenvector*eigenvalue < SMALL:" << nl;
        Info<< "First eigenvector:" << nl;
        for (direction i = 0; i < 2; ++i)
        {
            Vector2D<complex> evec;
            switch (i)
            {
                case 0: evec = EVecs.x(); break;
                case 1: evec = EVecs.y(); break;
            }
            Vector2D<complex> leftSide
            (
                T.xx()*evec.x() + T.xy()*evec.y(),
                T.yx()*evec.x() + T.yy()*evec.y()
            );
            Vector2D<complex> rightSide
            (
                evec.x()*EVals[i],
                evec.y()*EVals[i]
            );
            leftSide -= rightSide;
            isEqual(mag(leftSide), 0, verbose);
        }
    }

    Info<< nl;
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Random rndGen(1234);
    const label numberOfTests = 1000;

    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "### " << numberOfTests << " Symmetric Tensor2D<scalar>:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            tensor2D testTensor(randomTensor2D(rndGen));

            // Ensure testTensor is symmetric
            testTensor.xy() = testTensor.yx();

            Info<< "## testTensor = " << nl << testTensor << nl;

            verification_eigendecomposition_tensor2D(testTensor);
        }

        horizontalLine();
    }
    #endif

    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "### " << numberOfTests << " Nonsymm Tensor2D<scalar>:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            tensor2D testTensor(randomTensor2D(rndGen));

            Info<< "## testTensor = " << nl << testTensor << nl;

            verification_eigendecomposition_tensor2D(testTensor);
        }

        horizontalLine();
    }
    #endif

    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "### Repeating real eigenvalue examples:" << nl;
        {
            const tensor2D testTensor
            (
                2, 0,
                0, 2
            );

            verification_eigendecomposition_tensor2D(testTensor);
        }
        {
            // Defective matrix
            const tensor2D testTensor
            (
                3, 1,
                0, 3
            );

            verification_eigendecomposition_tensor2D(testTensor);
        }

        horizontalLine();
    }
    #endif


    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
