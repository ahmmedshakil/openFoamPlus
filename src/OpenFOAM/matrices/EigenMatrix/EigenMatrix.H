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

Class
    Foam::EigenMatrix

Description
    EigenMatrix (i.e. eigen decomposition or spectral decomposition) decomposes
    a diagonalisable nonsymmetric scalar square matrix into its canonical form,
    whereby the matrix is represented in terms of its eigenvalues and
    eigenvectors.
    
    The eigenvalue equation (i.e. eigenvalue problem) is:

    \verbatim
        A v = \lambda v
    \endverbatim

    where
     \c A is a diagonalisable square matrix of dimension m-by-m,
     \c v is a (non-zero) vector of dimension m (eigenvector),
     \c \lambda is a scalar corresponding to \c v (eigenvalue).

    Reference:
    \verbatim
        Eigen decomposition:
            Golub G. H., & Van Loan C. F. (2013).
            Matrix computations. 4th Edition.
            Baltimore: Johns Hopkins University Press.
            ISBN:9781421407944

            Bindel, D. (2012).
            Double-shift QR steps.
            Ithaca: Cornell University.
            bit.ly/2GFujbP

            Ford, W. (2014).
            Numerical linear algebra with applications: Using MATLAB.
            London: Elsevier/Academic Press.
            doi:10.1016/C2011-0-07533-6

            Watkins, D. S. (2002).
            Fundamentals of matrix computations.
            New York: John Wiley & Sons, Inc..
            doi:10.1002/0471249718

        Eigen decomposition properties:
            mathworld.wolfram.com/EigenDecomposition.html (Retrieved:15-06-19)
            mathworld.wolfram.com/EigenDecompositionTheorem.html 
            (Retrieved:15-06-19)

        Estimation for the maximum number of iterations in eigen decomposition:
            Kressner, D. (2005).
            Numerical methods for general and structured eigenvalue problems.
            Heidelberg: Springer Berlin Heidelberg.
            doi:10.1007/3-540-28502-4
    \endverbatim

Usage
    Input types:
    - Input \c A is a \c SquareMatrix<scalar>

    Output types:
    - Output \c EVals (eigenvalues) is a \c List<List<scalar>>
    - Output \c EVecs (eigenvectors) is a \c SquareMatrix<complex>

    Eigenvalues of \c EigenMatrix can be stored in the following ways:
    - storeMethod::IN_PLACE:      replaces input matrix content with EVal matrix
    - storeMethod::OUT_OF_PLACE:  creates new object of EVal matrix

    Computation of \c EVecs is optional by true/false Boolean or as follows:
     - computationEVecs::FALSE:   switches off the computation of \c EVecs
     - computationEVecs::TRUE:    switches on the computation of \c EVecs

Notes
    - \c EigenMatrix computes eigenvalues by the double-shift Francis algorithm
    using Householder reflectors in scalar arithmetic to avoid more expensive
    complex arithmetic. Yet output eigenvalues/eigenvectors can be complex.
    - \c EigenMatrix computes eigenvectors by the inverse power iteration.
    - \c EigenMatrix computes eigenvalues by default, and optionally all 
    or individual eigenvectors.
    - Any complex eigenvalue has its complex conjugate pair as an eigenvalue.
    - The name "Francis's algorithm" was used herein in line with Watkins (2002)
    since no explicit QR decompositions are performed in \c EigenMatrix.

See also
    Test-EigenMatrix.C

SourceFiles
    EigenMatrix.C
    EigenMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef EigenMatrix_H
#define EigenMatrix_H

#include "HessenbergMatrix.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class EigenMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class EigenMatrix:
    public HessenbergMatrix<MatrixType>
{

public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<scalar> SMatrix;
    typedef SymmetricSquareMatrix<scalar> SymMatrix;    // placeholder
    typedef RectangularMatrix<scalar> RMatrix;
    typedef SquareMatrix<complex> SCMatrix;
    typedef RectangularMatrix<complex> RCMatrix;

    static_assert
    (
        std::is_same<cmptType, scalar>::value,
        "EigenMatrix operates only in scalar arithmetic."
    );

    //- Options where to store EVals
    enum storeMethod : bool
    {
        IN_PLACE = false,    //!< replaces input matrix content with EVal matrix
        OUT_OF_PLACE = true  //!< creates new object of EVal matrix
    };

    //- Options for the computation of EVecs
    enum computationEVecs : bool
    {
        FALSE = false,      //!< switches off the computation of \c EVecs
        TRUE = true         //!< switches on the computation of \c EVecs
    };


private:

    // Private Data

        //- Selected option where to store EVals
        const storeMethod store_;

        //- Selected option for the computation of EVecs
        const computationEVecs computeEVecs_;

        //- Upper Hessenberg matrix that evolves to the eigenvalue matrix
        //  Diagonal contains either real EVals or real part of complex EVals
        //  Off-diagonal contains zeroes or imaginary part of complex EVals
        MatrixType& H_;

        //- Eigenvalues
        //  List[0] = Real EVals or real part of complex EVals
        //  List[1] = Zeroes or imaginary part of complex EVals
        List<List<scalar>> EVals_;

        //- Eigenvector matrix
        //  Each column is an EVec corresponds to an EVal
        SCMatrix EVecs_;


    // Private Member Functions

        //- Compute EVals by the Francis algorithm
        //  Left and right EVals are equivalent
        void eigenvalues
        (
            MatrixType& A,
            MatrixType& A0,
            const bool symmetric,
            const label maxIter,
            const scalar tol
        );

        //- Compute (in-place) implicit Francis double-shift
        //- stepping and bulge chasing
        void implicitStep
        (
            SMatrix& H
        );

        //- Return a Francis double-shift pair
        //  The shift pair is estimated from EVals of 2-by-2 trailing
        //  submatrices, and is used to improve the rate of convergence of
        //  iterative EVal computations
        Pair<scalar> FrancisDoubleShift
        (
            const scalar trace,
            const scalar determinant,
            const scalar A11
        ) const;

        //- Return EVals of a given 2-by-2 square matrix
        //  Return Pair(Eigenvalue1, Eigenvalue2) if real EVals
        //  Return Pair(Eigenvalue1Real, Eigenvalue1Imag) if complex EVals
        Pair<scalar> EVals2by2
        (
            const scalar trace,
            const scalar determinant,
            bool& isRealRoot
        ) const;

        //- Unpack real and imaginary parts of
        //- EVal matrix into a List of EVals
        void unpackEigenvalues
        (
            const MatrixType& A,
            const bool symmetric
        );

        //- Compute right EVecs corresponding to
        //- EVals by the inverse power iteration
        //  Left and right EVecs are different
        void eigenvectors
        (
            MatrixType& A0,
            const bool symmetric,
            const label maxIter,
            const scalar tol
        );


public:

    // Constructors

        //- Construct null
        EigenMatrix();

        //- Construct EigenMatrix without performing the decomposition
        EigenMatrix
        (
            const storeMethod store,
            const computationEVecs computeEVecs = computationEVecs::FALSE
        );

        //- Construct EigenMatrix and perform the eigen decomposition
        EigenMatrix
        (
            MatrixType& A,
            const storeMethod store,
            const computationEVecs computeEVecs = computationEVecs::FALSE,
            const label maxIterEVal = 2000,
            const scalar tolEVal = 1e-8,
            const label maxIterEVec = 10,
            const scalar tolEVec = 1e-8
        );


    // Member Functions

        // Information

        //- Return (complex) EVals
        inline const List<List<scalar>>& EVals() const;

        //- Return EVecs matrix
        inline const SCMatrix& EVecs() const;


        // Algorithm

        //- Compute eigen decomposition according to the chosen settings
        //  (Golub and Van Loan, (2013), Algorithm 7.5.1 and 7.5.2)
        void decompose
        (
            MatrixType& A,
            label maxIterEVal = -1,
            scalar tolEVal = -1,
            label maxIterEVec = -1,
            scalar tolEVec = -1
        );

        //- Compute the right EVec corresponding to the given EVal
        RMatrix eigenvector
        (
            MatrixType& A0,
            const scalar eVal,
            const label maxIter = 10,
            const scalar tol = 1e-8,
            const scalar perturb = 0
        );

        //- Compute the right EVec corresponding to the given EVal
        RCMatrix eigenvector
        (
            const MatrixType& A0,
            const complex& eVal,
            const label maxIter = 10,
            const scalar tol = 1e-8,
            const scalar perturb = 0
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EigenMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EigenMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
