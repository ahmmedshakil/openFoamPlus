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
#include "complex.H"
#include "Random.H"
#include "HessenbergMatrix.H"
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


// Create random scalar-type matrix
template<class MatrixType>
typename std::enable_if
<
    !std::is_same<complex, typename MatrixType::cmptType>::value,
    MatrixType
>::type makeRandomMatrix
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


// Create random complex-type matrix
template<class MatrixType>
typename std::enable_if
<
    std::is_same<complex, typename MatrixType::cmptType>::value,
    MatrixType
>::type makeRandomMatrix
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
        [&]
        {
            return complex
            (
                rndGen.GaussNormal<scalar>(),
                rndGen.GaussNormal<scalar>()
            );
        }
    );

    return mat;
}


// Print OpenFOAM matrix in NumPy format
template<class MatrixType>
void InfoNumPyFormat
(
    const MatrixType& mat,
    const word objName
)
{
    Info<< objName << ": " << mat.m() << "x" << mat.n() << nl;
    for (label m = 0; m < mat.m(); ++m)
    {
        Info<< "[";
        for (label n = 0; n < mat.n(); ++n)
        {
            if (n == mat.n() - 1)
            {
                Info<< mat(m,n);
            }
            else
            {
                Info<< mat(m,n) << ",";
            }
        }
        if (m == mat.m() - 1)
        {
            Info << "]" << nl;
        }
        else
        {
            Info << "]," << nl;
        }
    }
}


// Returns true if two scalars are equal within a given tolerance, and v.v.
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


// Checks if H matrix is an upper Hessenberg matrix
template<class Type>
void cross_check_HessenbergMatrix
(
    const SquareMatrix<Type>& Aorig,
    const SquareMatrix<Type>& H
)
{
    if (Aorig.symmetric())
    {
        Info<< "Tridiagonal Hessenberg = " << H.tridiagonal() << nl;
    }

    InfoNumPyFormat(Aorig, "Aorig");
    InfoNumPyFormat(H, "H");

    // mathworld.wolfram.com/HessenbergMatrix.html (Retrieved:13-06-19)
    Info<< nl << "# H(i, j) = 0 for i > j + 1:" << nl;
    for (label i = 0; i < H.m(); ++i)
    {
        for (label j = 0; j < H.n(); ++j)
        {
            if (j + 1 < i)
            {
                isEqual(mag(H(i, j)), 0.0, verbose);
            }
        }
    }
}


// Checks if Q matrix is a unitary matrix
// Checks if the given matrix can be reconstructed by H and Q matrices
template<class Type>
void cross_check_HessenbergMatrix
(
    const SquareMatrix<Type>& Aorig,
    const SquareMatrix<Type>& H,
    const SquareMatrix<Type>& Q
)
{
    InfoNumPyFormat(Q, "Q");

    // mathworld.wolfram.com/HessenbergDecomposition.html (Retrieved:13-06-19)
    {
        Info<< nl << "# Original Input A = Q*H*Q.T():" << nl;
        const SquareMatrix<Type> Areconstruct(Q*(H^Q));
        equal(Aorig, Areconstruct, verbose);
    }

    // mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:16-06-19)
    {
        Info<< nl << "# Q.T() = Q^-1:" << nl;
        const SquareMatrix<Type> Qpinv(pinv(Q));
        equal(Q.T(), Qpinv, verbose);
    }

    // mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:16-06-19)
    {
        Info<< nl << "# Q.T()*Q = Q*Q.T() = I:" << nl;
        const SquareMatrix<Type> QTQ(Q&Q);
        const SquareMatrix<Type> QQT(Q^Q);
        const SquareMatrix<Type> IMatrix(Q.m(), I);
        equal(QTQ, IMatrix, verbose);
        equal(QQT, IMatrix, verbose);
        equal(QTQ, QQT, verbose);
    }

    cross_check_HessenbergMatrix(Aorig, H);
}


// Checks each constructor of HessenbergMatrix type
template<class Type>
void verification_HessenbergMatrix
(
    SquareMatrix<Type>& A
)
{
    typedef SquareMatrix<Type> SMatrix;

    // Create a copy of matrix A
    const SMatrix Aorig(A);

    // HessenbergMatrix Constructors
    #if (0 | RUNALL)
    {
        Info<< "# Null constructor" << nl;
        HessenbergMatrix<SMatrix> HMNull();
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE " << nl;
        HessenbergMatrix<SMatrix> HM
        (
            HessenbergMatrix<SMatrix>::storeMethod::IN_PLACE
        );
        SMatrix A0(A);
        HM.decompose(A0);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(H, "H");
        isEmpty(Q, "Q");
        cross_check_HessenbergMatrix(Aorig, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# IN_PLACE | computeQ = true" << nl;
        HessenbergMatrix<SMatrix> HM
        (
            HessenbergMatrix<SMatrix>::storeMethod::IN_PLACE,
            HessenbergMatrix<SMatrix>::computationQ::TRUE
        );
        SMatrix A0(A);
        HM.decompose(A0);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(H, "H");
        cross_check_HessenbergMatrix(Aorig, A0, Q);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | IN_PLACE " << nl;
        SMatrix A0(A);
        HessenbergMatrix<SMatrix> HM
        (
            A0,
            HessenbergMatrix<SMatrix>::storeMethod::IN_PLACE
        );
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(H, "H");
        isEmpty(Q, "Q");
        cross_check_HessenbergMatrix(Aorig, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | IN_PLACE | computeQ = true" << nl;
        SMatrix A0(A);
        HessenbergMatrix<SMatrix> HM
        (
            A0,
            HessenbergMatrix<SMatrix>::storeMethod::IN_PLACE,
            HessenbergMatrix<SMatrix>::computationQ::TRUE
        );
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(H, "H");
        cross_check_HessenbergMatrix(Aorig, A0, Q);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE " << nl;
        HessenbergMatrix<SMatrix> HM
        (
            HessenbergMatrix<SMatrix>::storeMethod::OUT_OF_PLACE
        );
        SMatrix A0(A);
        HM.decompose(A0);
        equal(A0, A, verbose);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(Q, "Q");
        cross_check_HessenbergMatrix(Aorig, H);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# OUT_OF_PLACE | computeQ = true" << nl;
        HessenbergMatrix<SMatrix> HM
        (
            HessenbergMatrix<SMatrix>::storeMethod::OUT_OF_PLACE,
            HessenbergMatrix<SMatrix>::computationQ::TRUE
        );
        SMatrix A0(A);
        HM.decompose(A0);
        equal(A0, A, verbose);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        cross_check_HessenbergMatrix(Aorig, H, Q);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | OUT_OF_PLACE " << nl;
        SMatrix A0(A);
        HessenbergMatrix<SMatrix> HM
        (
            A0,
            HessenbergMatrix<SMatrix>::storeMethod::OUT_OF_PLACE
        );
        equal(A0, A, verbose);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        isEmpty(Q, "Q");
        cross_check_HessenbergMatrix(Aorig, H);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | OUT_OF_PLACE | computeQ = true" << nl;
        SMatrix A0(A);
        HessenbergMatrix<SMatrix> HM
        (
            A0,
            HessenbergMatrix<SMatrix>::storeMethod::OUT_OF_PLACE,
            HessenbergMatrix<SMatrix>::computationQ::TRUE
        );
        equal(A0, A, verbose);
        const SMatrix& H(HM.H());
        const SMatrix& Q(HM.Q());
        cross_check_HessenbergMatrix(Aorig, H, Q);
    }
    #endif
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    typedef SquareMatrix<scalar> SMatrix;
    typedef SquareMatrix<complex> SCMatrix;

    Info<< setprecision(15);
    Random rndGen(12345);
    label numberOfTests = 100;

    Info<< "### Hessenberg decomposition:" << nl << nl;

    #if (0 | RUNALL)
    // SquareMatrix<scalar>
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<scalar> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(1, 50);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));
            verification_HessenbergMatrix(A);
        }

        horizontalLine();
    }
    #endif

    // SquareMatrix<complex>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<complex> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position(1, 50);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SCMatrix A(makeRandomMatrix<SCMatrix>({mRows, mRows}, rndGen));
            verification_HessenbergMatrix(A);
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
            const label mRows = rndGen.position(1, 50);
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
            verification_HessenbergMatrix(A);
        }

        horizontalLine();
    }
    #endif


    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
