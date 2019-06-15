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

#include "HessenbergMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
void Foam::HessenbergMatrix<MatrixType>::hessenberg
(
    MatrixType& A
)
{
    #ifdef FULLDEBUG
    if (A.empty())
    {
        FatalErrorInFunction
            << "Empty input Matrix" << abort(FatalError);
    }
    #endif

    const label mRows = A.m();

    // Allocate resources for Q_, if need be
    if (computeQ_)
    {
        Q_ = SMatrix(mRows, I);
    }

    // No action for 2-by-2 and 1-by-1 matrices
    // For 3-by-3 matrices, check if the leftmost-bottom elem is already zero
    if (mRows < 3)
    {
        return;
    }
    else if (mRows == 3 && (mag(A(2, 0))) < SMALL)
    {
        A(2, 0) = cmptType(0);
        return;
    }

    for (label k = 0; k < mRows - 2; ++k)
    {
        const label k1 = k + 1;

        // Note that Householder reflector selection is different for QRMatrix
        const RMatrix reflector
        (
            QRMatrix<MatrixType>::householderReflector(A.subColumn(k, k1))
        );

        // Unitary similarity matrix
        if (computeQ_)
        {
            QRMatrix<MatrixType>::applyRightReflector(Q_, reflector, k1);
        }

        // u*A; Operate on the entire matrix except k^th row and its above
        QRMatrix<MatrixType>::applyLeftReflector(A, reflector, k, k1);

        // (u*A)*u; Operate on the entire matrix except k^th column and its left
        QRMatrix<MatrixType>::applyRightReflector(A, reflector, k1);
    }

    A.round();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MatrixType>
void Foam::HessenbergMatrix<MatrixType>::decompose
(
    MatrixType& A
)
{
    switch (store_)
    {
        case storeMethod::IN_PLACE:
            hessenberg(A);
            break;

        case storeMethod::OUT_OF_PLACE:
            H_ = A;
            hessenberg(H_);
            break;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::HessenbergMatrix<MatrixType>::HessenbergMatrix()
{}


template<class MatrixType>
Foam::HessenbergMatrix<MatrixType>::HessenbergMatrix
(
    const storeMethod store,
    const computationQ computeQ
)
:
    QRMatrix<MatrixType>(),
    store_(store),
    computeQ_(computeQ),
    H_(),
    Q_()
{}


template<class MatrixType>
Foam::HessenbergMatrix<MatrixType>::HessenbergMatrix
(
    MatrixType& A,
    const storeMethod store,
    const computationQ computeQ
)
:
    QRMatrix<MatrixType>(),
    store_(store),
    computeQ_(computeQ),
    H_(),
    Q_()
{
    decompose(A);
}


// ************************************************************************* //
