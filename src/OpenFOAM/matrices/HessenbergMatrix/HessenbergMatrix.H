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
    Foam::HessenbergMatrix

Description
    HessenbergMatrix (i.e. Hessenberg decomposition) decomposes a square matrix
    \c A into an upper Hessenberg matrix \c H, and a unitary similarity matrix
    \c Q as follows:

    \verbatim
        A = Q H Q^(T)
    \endverbatim

    where \c T is the (conjugate) transpose.

    \c HessenbergMatrix can be computed in a fixed number of steps, and helps
    to accelerate by-nature iterative eigen decomposition algorithms since
    \c HessenbergMatrix and original matrix will have the same eigenvalues.

    Reference:
    \verbatim
        Hessenberg decomposition:
            Golub G. H., & Van Loan C. F. (2013).
            Matrix computations. 4th Edition.
            Baltimore: Johns Hopkins University Press.
            ISBN:9781421407944

        Operands of Hessenberg decomposition:
            mathworld.wolfram.com/HessenbergMatrix.html (Retrieved:15-06-19)
            mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:13-06-19)
    \endverbatim

Usage
    Input types:
     - Input \c A is a \c SquareMatrix<Type>

    Output types:
     - Output \c H (upper Hessenberg matrix) is a \c SquareMatrix<Type>
     - Output \c Q (unitary similarity matrix) is a \c SquareMatrix<Type>

    \c H of \c HessenbergMatrix can be stored in the following ways:
     - storeMethod::IN_PLACE:        replaces input matrix content with \c H
     - storeMethod::OUT_OF_PLACE:    creates new object of \c H

    Computation of \c Q is optional by true/false Boolean or as follows:
     - computationQ::FALSE:          switches off the computation of Q
     - computationQ::TRUE:           switches on the computation of Q 

Notes
    - Hessenberg decomposition is not unique (ibid., p. 381).
    - (Conjugate) transpose of Hessenberg matrix is also a Hessenberg matrix.
    - A Hessenberg-decomposed symmetric square matrix is also symmetric.
    - Hessenberg decomposition reduces a Hermitian matrix to a tridiagonal form.

See also
    Test-HessenbergMatrix.C

SourceFiles
    HessenbergMatrix.C
    HessenbergMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef HessenbergMatrix_H
#define HessenbergMatrix_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "QRMatrix.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HessenbergMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class HessenbergMatrix:
    public QRMatrix<MatrixType>
{

public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> SMatrix;
    typedef SymmetricSquareMatrix<cmptType> SymMatrix;      // placeholder
    typedef RectangularMatrix<cmptType> RMatrix;

    static_assert
    (
        std::is_same<MatrixType, SMatrix>::value
     || std::is_same<MatrixType, SymMatrix>::value,
        "Hessenberg decomposition is not defined for non-square matrices."
    );

    static_assert
    (
        !std::is_same<MatrixType, SymMatrix>::value,
        "Symmetric square matrices are not supported, yet."
    );

    //- Options where to store H
    enum storeMethod : bool
    {
        IN_PLACE = false,      //!< replaces input matrix content with \c H
        OUT_OF_PLACE = true    //!< creates new object of \c H
    };

    //- Options for the computation of Q
    enum computationQ : bool
    {
        FALSE = false,      //!< switches off the computation of Q
        TRUE = true         //!< switches on the computation of Q
    };


private:

    // Private Member Functions

        //- (In-place) compute (upper) Hessenberg decomposition
        //  (Golub and Van Loan, (2013), Algorithm 7.4.2.)
        void hessenberg
        (
            MatrixType& A
        );


protected:

    // Protected Data

        //- Selected option where to store H
        const storeMethod store_;

        //- Selected option for the computation of Q
        const computationQ computeQ_;

        //- Upper Hessenberg matrix
        MatrixType H_;

        //- Unitary similarity matrix
        MatrixType Q_;


public:

    // Constructors

        //- Construct null
        HessenbergMatrix();

        //- Construct HessenbergMatrix without performing the decomposition
        HessenbergMatrix
        (
            const storeMethod store,
            const computationQ computeQ = computationQ::FALSE
        );

       //- Construct HessenbergMatrix and perform the Hessenberg decomposition
       HessenbergMatrix
       (
           MatrixType& A,
           const storeMethod store,
           const computationQ computeQ = computationQ::FALSE
       );


    // Member Functions

        // Information

        //- Return the upper Hessenberg matrix
        //  Return an empty MatrixType if storeMethod::IN_PLACE
        inline const MatrixType& H() const;

        //- Return the unitary similarity matrix
        //  Return an empty MatrixType if computeQ_=false
        inline const MatrixType& Q() const;


        // Algorithm

        //- Compute Hessenberg decomposition according to constructor settings
        void decompose
        (
            MatrixType& A
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HessenbergMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HessenbergMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
