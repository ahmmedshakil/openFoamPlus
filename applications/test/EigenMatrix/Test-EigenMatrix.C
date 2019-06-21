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

#include "MatrixTools.H"
#include "scalarMatrices.H"
#include "Random.H"
#include "EigenMatrix.H"
#include "IOmanip.H"

using namespace Foam;
using namespace Foam::MatrixTools;

#define equal MatrixTools::equal
#define RUNALL true
const bool verbose = true;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void horizontalLine()
{
    Info<< "+---------+---------+---------+---------+---------+" << nl;
}


// Copy values into matrix
template<class Form, class Type>
void assignMatrix
(
    Matrix<Form, Type>& mat,
    std::initializer_list<typename Matrix<Form, Type>::cmptType> list
)
{
    const label nargs = list.size();

    if (nargs != mat.size())
    {
        FatalErrorInFunction
            << "Mismatch in matrix dimension ("
            << mat.m() << ", "
            << mat.n() << ") and number of args (" << nargs << ')' << nl
            << exit(FatalError);
     }

    std::copy(list.begin(), list.end(), mat.begin());
}


// Create random scalar-type matrix
template<class MatrixType>
MatrixType makeRandomMatrix
(
    const labelPair& dims,
    Random& rndGen
)
{
    MatrixType mat(dims);

    std::generate
    (
        mat.begin(),
        mat.end(),
        [&]{return rndGen.GaussNormal<typename MatrixType::cmptType>();}
    );

    return mat;
}


// Print OpenFOAM matrix in NumPy format
template<class MatrixType>
void InfoNumPyFormat
(
    const MatrixType& mat,
    const bool matlab = false
)
{
    for (label m = 0; m < mat.m(); ++m)
    {
        if (matlab == false || m == 0)
        {
            Info<< "[";
        }
        for (label n = 0; n < mat.n(); ++n)
        {
            if (n == mat.n() - 1)
            {
                Info<< mat(m,n);
            }
            else
            {
                if (matlab == true)
                {
                    Info<< mat(m,n) << " ";
                }
                else
                {
                    Info<< mat(m,n) << ",";
                }
            }
        }
        if (m == mat.m() - 1)
        {
            Info << "]" << nl;
        }
        else
        {
            if (matlab == true)
            {
                Info << ";" << nl;
            }
            else
            {
                Info << "]," << nl;
            }
        }
    }
}


// Returns true if two scalars are equal within a given tolerance, and v.v.
bool isEqual
(
    const scalar a,
    const scalar b,
    const bool verbose = false,
    const scalar relTol = 1e-3,
    const scalar absTol = 1e-3
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


// Prints true if a given matrix is empty, and v.v.
template<class MatrixType>
void isEmpty
(
    const MatrixType& A,
    const word objName
)
{
    Info<< "Empty " << objName << " = ";
    if (A.empty())
    {
        Info<< "true" << nl;
    }
    else
    {
        Info<< "false" << nl;
    }
}


// Checks if eigenvalues satify eigen decomposition relations
template<class MatrixType>
void cross_check_eigenvalues
(
    const MatrixType& Aorig,
    const List<scalar>& EValsReal,
    const List<scalar>& EValsImag
)
{
    // w.wiki/4zs (Retrieved: 16-06-19) # Item-1
    #if (1 | RUNALL)
    {
        Info<< nl << "# Sum(Eigenvalues) = trace(A):" << nl;
        const scalar trace = Aorig.trace();
        // Imaginary part of complex conjugates cancel each other
        const scalar EValsSum = sum(EValsReal);
        isEqual(EValsSum, trace, verbose);
    }
    #endif

    // w.wiki/4zs (Retrieved: 16-06-19) # Item-2
    #if (1 | RUNALL)
    {
        Info<< nl << "# Prod(Eigenvalues) = det(A):"
            << " (Note that the det computation may fail.)" << nl;
        const scalar determinant = mag(det(Aorig));
        scalar EValsProd = 1.0;

        if (EValsImag.empty())
        {
            for (label i = 0; i < EValsReal.size(); ++i)
            {
                EValsProd *= Foam::sqrt(sqr(EValsReal[i]));
            }
        }
        else
        {
            for (label i = 0; i < EValsReal.size(); ++i)
            {
                EValsProd *= Foam::sqrt(sqr(EValsReal[i]) + sqr(EValsImag[i]));
            }
        }
        isEqual(EValsProd, determinant, verbose);
    }
    #endif
}


// Checks if eigenvectors satify eigen decomposition relations
template<class MatrixType>
void cross_check_eigenvectors
(
    const MatrixType& Aorig,
    const List<scalar>& EValsReal,
    const List<scalar>& EValsImag,
    const SquareMatrix<complex>& EVecs
)
{
    #if (1 | RUNALL)
    {
        Info<< nl << "# (A*EVec - EVal*EVec) < SMALL" << nl;
        Info<< nl << "# (A - Eigenvalues*I)*Eigenvectors = 0" << nl;

        SquareMatrix<complex> A(Aorig.m());
        auto convertToComplex = [&](const scalar& val) { return complex(val); };
        std::transform
        (
            Aorig.cbegin(),
            Aorig.cend(),
            A.begin(),
            convertToComplex
        );

        for (label i = 0; i < Aorig.m(); ++i)
        {
            const RectangularMatrix<complex>& EVec(EVecs.subColumn(i));
            const complex EVal(EValsReal[i], EValsImag[i]);
            const RectangularMatrix<complex> leftSide(A*EVec);
            const RectangularMatrix<complex> rightSide(EVal*EVec);
            equal(leftSide, rightSide, verbose, 10, 1e-3, 1e-3);
        }
    }
    #endif
}


// Checks if eigenvalues and eigenvectors satify eigen decomposition relations
template<class MatrixType>
void cross_check_EigenMatrix
(
    const MatrixType& Aorig,
    EigenMatrix<MatrixType>& EM,
    const bool computeEVecs = false
)
{
    typedef SquareMatrix<complex> SCMatrix;

    Info<< "Aorig" << nl;
    const bool infoMatlabFormat = false;
    InfoNumPyFormat(Aorig, infoMatlabFormat);

    // Eigenvalue checks
    const List<List<scalar>>& EValsList(EM.EVals());
    Info<< "EVals:" << nl << EValsList << nl;
    cross_check_eigenvalues(Aorig, EValsList[0], EValsList[1]);

    // Eigenvector checks
    if (computeEVecs)
    {
        const SCMatrix& EVecs(EM.EVecs());
        Info<< "EVecs:" << nl << EVecs << nl;
        cross_check_eigenvectors(Aorig, EValsList[0], EValsList[1], EVecs);
    }
}


// Checks each constructor of EigenMatrix type
template<class MatrixType>
void verification_EigenMatrix
(
    MatrixType& A
)
{
    // Create a copy of matrix A
    const MatrixType Aorig(A);
    const scalar maxIter = 2000;
    const scalar tol = 1e-5;

    // EigenMatrix Constructors
    #if (0 | RUNALL)
    {
        Info<< "# Null constructor" << nl;
        EigenMatrix<MatrixType> EMNull();
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE " << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            EigenMatrix<MatrixType>::storeMethod::IN_PLACE
        );
        EM.decompose(A0, maxIter, tol);
        cross_check_EigenMatrix(Aorig, EM);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE | computeEVecs = false | maxIter | tol" << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            A0,
            EigenMatrix<MatrixType>::storeMethod::IN_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::FALSE,
            maxIter,
            tol
        );
        cross_check_EigenMatrix(Aorig, EM);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE " << nl;
        EigenMatrix<MatrixType> EM
        (
            EigenMatrix<MatrixType>::storeMethod::OUT_OF_PLACE
        );
        EM.decompose(A, maxIter, tol);
        cross_check_EigenMatrix(Aorig, EM);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE | computeEVecs = false | maxIter | tol" << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            A0,
            EigenMatrix<MatrixType>::storeMethod::OUT_OF_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::FALSE,
            maxIter,
            tol
        );
        cross_check_EigenMatrix(Aorig, EM);
    }
    #endif

    // Eigenvector computation
    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE | computeEVecs = true" << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            EigenMatrix<MatrixType>::storeMethod::IN_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::TRUE
        );
        EM.decompose(A0, maxIter, tol);
        cross_check_EigenMatrix(Aorig, EM, verbose);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE | computeEVecs = true | maxIter | tol" << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            A0,
            EigenMatrix<MatrixType>::storeMethod::IN_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::TRUE,
            maxIter,
            tol
        );
        cross_check_EigenMatrix(Aorig, EM, verbose);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE | computeEVecs = true" << nl;
        EigenMatrix<MatrixType> EM
        (
            EigenMatrix<MatrixType>::storeMethod::OUT_OF_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::TRUE
        );
        EM.decompose(A, maxIter, tol);
        cross_check_EigenMatrix(Aorig, EM, verbose);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE | computeEVecs = true | maxIter | tol" << nl;
        MatrixType A0(A);
        EigenMatrix<MatrixType> EM
        (
            A0,
            EigenMatrix<MatrixType>::storeMethod::OUT_OF_PLACE,
            EigenMatrix<MatrixType>::computationEVecs::TRUE,
            maxIter,
            tol
        );
        cross_check_EigenMatrix(Aorig, EM, verbose);
    }
    #endif
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    typedef SquareMatrix<scalar> SMatrix;
    Info<< setprecision(15);
    Random rndGen(1234);
    label numberOfTests = 20;

    Info<< "### Eigen decomposition:" << nl << nl;

    // SquareMatrix<scalar>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<scalar> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(1, 100);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));

            verification_EigenMatrix(A);
        }

        horizontalLine();
    }
    #endif

    // Symmetric SquareMatrix<scalar>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " Symm. SquareMatrix<scalar> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(1, 100);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));
            // Symmetrise, so that eigenvalues & eigenvectors are real
            for (label n = 0; n < A.n() - 1; ++n)
            {
                for (label m = A.m() - 1; m > n; --m)
                {
                    A(n, m) = A(m, n);
                }
            }

            verification_EigenMatrix(A);
        }

        horizontalLine();
    }
    #endif

    #if (0 | RUNALL)
    {
        // Rosser, J. B., Lanczos, C., Hestenes, M. R., & Karush, W. (1951).
        // Separation of close eigenvalues of a real symmetric matrix.
        // Journal of Research of the National Bureau of Standards, 47(4),
        // 291-297. doi: 10.6028/jres.047.037
        // 
        // 8x8 symmetric square matrix consisting of close real eigenvalues
        // ibid, p. 294
        // {
        //      1020.04901843, 1020.000, 1019.90195136, 1000.000,
        //      1000.000, 0.09804864072, 0.000, -1020.0490
        // }
        // Note that Prod(Eigenvalues) != determinant(A) for this matrix
        // via the LAPACK routine z/dgetrf

        horizontalLine();
        SMatrix A({8, 8}, Zero);
        assignMatrix
        (
            A,
            {
                611, 196, -192, 407, -8, -52, -49, 29,
                196, 899, 113, -192, -71, -43, -8, -44,
                -192, 113, 899, 196, 61, 49, 8, 52,
                407, -192, 196, 611, 8, 44, 59, -23,
                -8, -71, 61, 8, 411, -599, 208, 208,
                -52, -43, 49, 44, -599, 411, 208, 208,
                -49, -8, 8, 59, 208, 208, 99, -911,
                29, -44, 52, -23, 208, 208, -911, 99
            }
        );

        verification_EigenMatrix(A);

        horizontalLine();
    }
    #endif

    // Repeating real eigenvalues
    #if (0 | RUNALL)
    {
        {
            SMatrix A({3, 3}, Zero);
            assignMatrix
            (
                A,
                {
                    0, 1, 1,
                    1, 0, 1,
                    1, 1, 0
                }
            );
            verification_EigenMatrix(A);
        }
        {
            SMatrix A({3, 3}, Zero);
            assignMatrix
            (
                A,
                {
                    2, 0, 0,
                    0, 2, 0,
                    0, 0, 1
                }
            );
            verification_EigenMatrix(A);
        }
    }
    #endif


    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
