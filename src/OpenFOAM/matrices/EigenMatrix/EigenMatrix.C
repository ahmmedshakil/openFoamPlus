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

#include "EigenMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
void Foam::EigenMatrix<MatrixType>::eigenvalues
(
    MatrixType& A,
    MatrixType& A0,
    const bool symmetric,
    const label maxIter,
    const scalar tol
)
{
    if (computeEVecs_)
    {
        A0 = A;
    }

    label iter = 0;
    label nCols = A.n();

    while (2 < nCols)
    {
        // Deflate 1-by-1 block
        if (mag(A(nCols - 1, nCols - 2)) < tol)
        {
            A(nCols - 1, nCols - 2) = 0.0;
            --nCols;
        }
        // Deflate 2-by-2 block if matrix is nonsymmetric
        else if (!symmetric && mag(A(nCols - 2, nCols - 3)) < tol)
        {
            A(nCols - 2, nCols - 3) = 0.0;
            nCols -= 2;
        }
        else
        {
            SMatrix H0(A.subMatrix(0, 0, nCols, nCols));
            implicitStep(H0);
            A.subMatrix(0, 0, nCols, nCols) = H0;
        }

        ++iter;
        if (maxIter < iter)
        {
            WarningInFunction
                << "Francis algorithm did not converge within: "
                << iter << " iterations." << nl;
            break;
        }
    }

    // Last 2-by-2 subMatrix or Matrix itself
    if (nCols == 2)
    {
        bool isRealRoot = true;
        const Pair<scalar> eVals
        (
            EVals2by2
            (
                A(0, 0) + A(1, 1),                    // trace
                A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0),    // determinant
                isRealRoot
            )
        );

        if (isRealRoot)
        {
            A(0, 0) = eVals.first();
            A(1, 1) = eVals.second();
            A(0, 1) = Zero;
            A(1, 0) = Zero;
        }
        else
        {
            A(0, 0) = eVals.first();
            A(1, 1) = A(0, 0);
            A(0, 1) = eVals.second();
            A(1, 0) = -A(0, 1);
        }
    }

    A.round();

    unpackEigenvalues(A, symmetric);
}


template<class MatrixType>
void Foam::EigenMatrix<MatrixType>::implicitStep
(
    SMatrix& H
)
{
    const label mRows = H.m();
    const label nCols = H.n();

    RMatrix C(3, mRows, Zero);
    {
        // 2-by-2 submatrix at the bottom-right corner of H
        const label m0 = mRows - 2;
        const label m1 = mRows - 1;
        const Pair<scalar> shifts
        (
            FrancisDoubleShift
            (
                H(m0, m0) + H(m1, m1),                        // trace
                H(m0, m0)*H(m1, m1) - H(m0, m1)*H(m1, m0),    // determinant
                H(m1, m1)
            )
        );

        // C = H.subMatrix(0, 0, 3, 2)*H.subColumn(0, 0, 2); 3-by-2*2-by-mRows
        for (label i = 0; i < C.m(); ++i)
        {
            for (label k = 0; k < 2; ++k)
            {
                for (label j = 0; j < C.n(); ++j)
                {
                    C(i, j) += H(i, k)*H(k, j);
                }
            }
        }

        // C.subColumn(0, 0, 2) += shifts.first()*H1; 2-by-mRows
        for (label i = 0; i < 2; ++i)
        {
            for (label j = 0; j < C.n(); ++j)
            {
                C(i, j) += shifts.first()*H(i, j);
            }
        }
        C(0, 0) += shifts.second();
        C(0, 0) += C(0, 0)/mag(C(0, 0))*C.columnNorm(0);
        C /= cmptType(C.columnNorm(0));
    }

    // Application of Householder reflectors
    {
        QRMatrix<MatrixType>::applyLeftReflector(H, C);

        QRMatrix<MatrixType>::applyRightReflector(H, C);
    }

    // Application of 'Bulge Chasing'
    for(label j = 0; j < nCols - 2; ++j)
    {
        const label k = min(j + 3, nCols - 1);

        const RMatrix u
        (
            QRMatrix<MatrixType>::householderReflector
            (
                H.subColumn(j, j + 1, k - (j + 1) + 1)
            )
        );

        QRMatrix<MatrixType>::applyLeftReflector(H, u, 0, j + 1);

        QRMatrix<MatrixType>::applyRightReflector(H, u, j + 1);

        H(k, j) = 0.0;
    }
}


template<class MatrixType>
Foam::Pair<Foam::scalar> Foam::EigenMatrix<MatrixType>::FrancisDoubleShift
(
    const scalar trace,
    const scalar determinant,
    const scalar A11
) const
{
    bool isRealRoot = true;
    Pair<scalar> eVals
    (
        EVals2by2
        (
            trace,
            determinant,
            isRealRoot
        )
    );

    // Eigenvalues are real
    if (isRealRoot)
    {
        if (mag(eVals.first() - A11) < mag(eVals.second() - A11))
        {
            eVals.second() = eVals.first();
        }
        else
        {
            eVals.first() = eVals.second();
        }

        return Pair<scalar>
        (
           -eVals.first() - eVals.second(),        // shift1
            eVals.first()*eVals.second()            // shift2
        );
    }
    // Eigenvalues are complex, (Watkins, (2002), p. 388)
    else
    {
        return Pair<scalar>
        (
            -trace,
            determinant
        );
    }
}


template<class MatrixType>
Foam::Pair<Foam::scalar> Foam::EigenMatrix<MatrixType>::EVals2by2
(
    const scalar trace,
    const scalar determinant,
    bool& isRealRoot
) const
{
    // Square-distance between two eigenvalues
    const scalar gapSqr = sqr(trace) - 4.0*determinant;

    // (Ford, (2014))
    // Eigenvalues are real
    if (0 <= gapSqr)
    {
        isRealRoot = true;
        const scalar firstRoot = 0.5*(trace - sign(-trace)*Foam::sqrt(gapSqr));

        #ifdef FULLDEBUG
            if (mag(firstRoot) < SMALL)
            FatalErrorInFunction
                << "SMALL-valued root is detected."
                << abort(FatalError);
        #endif

        return Pair<scalar>
        (
            firstRoot,
            determinant/firstRoot
        );
    }
    //Eigenvalues are complex conjugate pairs, return single pair member
    else
    {
        isRealRoot = false;
        return Pair<scalar>
        (
            0.5*trace,
            0.5*Foam::sqrt(mag(gapSqr))
        );
    }
}


template<class MatrixType>
void Foam::EigenMatrix<MatrixType>::unpackEigenvalues
(
    const MatrixType& A,
    const bool symmetric
)
{
    label i = 0;
    EVals_[0] = A.diag();
    EVals_[1] = List<scalar>(EVals_[0].size(), Zero);

    if (symmetric == false)
    {
        const label rowOffset = 1;
        List<cmptType> lowerSubDiag(min(A.m() - rowOffset, A.n()));

        label k = 0;
        for (cmptType& val : lowerSubDiag)
        {
            val = A(k + rowOffset, k);
            ++k;
        }

        while (i < lowerSubDiag.size())
        {
            if (SMALL < mag(lowerSubDiag[i]))
            {
                bool isRealRoot = true;
                const Pair<scalar> eVals
                (
                    EVals2by2
                    (
                        A(i, i) + A(i + 1, i + 1),
                        A(i, i)*A(i + 1, i + 1) - A(i, i + 1)*A(i + 1, i),
                        isRealRoot
                    )
                );

                EVals_[0][i] = eVals.first();
                EVals_[0][i+1] = eVals.first();
                EVals_[1][i] = eVals.second();
                EVals_[1][i+1] = -eVals.second();
                i += 2;
            }
            else
            {
                ++i;
            }
        }
    }
}


template<class MatrixType>
void Foam::EigenMatrix<MatrixType>::eigenvectors
(
    MatrixType& A0,
    const bool symmetric,
    const label maxIter,
    const scalar tol
)
{
    const label mRows = A0.m();
    EVecs_.resize(mRows);

    const MatrixType& Q = HessenbergMatrix<MatrixType>::Q_;

    if (symmetric)
    {
        for (label k = 0; k < mRows; ++k)
        {
            RMatrix eVec(eigenvector(A0, EVals_[0][k], maxIter, tol));

            // (Ford, (2014), Section 18.10)
            // Check NaNs and Infs
            bool hasNan = false;
            label p = 0;
            while (p < eVec.m() && !hasNan)
            {
                if (std::isnan(eVec(p, 0)) || std::isinf(eVec(p, 0)))
                {
                    hasNan = true;
                }
                ++p;
            }

            // In case of NaNs or Infs, perturb EVal and redo EVec computation
            if (hasNan)
            {
                Random rndGen(1234);
                const scalar perturb
                (
                    (2 + mag(rndGen.GaussNormal<scalar>()))*SMALL
                );

                eVec = eigenvector(A0, EVals_[0][k], maxIter, tol, perturb);
            }

            //EVecs_.subColumn(k) = Q*eVec;
            for (label i = 0; i < EVecs_.m(); ++i)
            {
                for (label q = 0; q < eVec.m(); ++q)
                {
                    EVecs_(i, k).Re() += Q(i, q)*eVec(q, 0);
                }
            }
        }
    }
    else
    {
        for (label k = 0; k < mRows; ++k)
        {
            const complex eVal(EVals_[0][k], EVals_[1][k]);
            RCMatrix eVec(eigenvector(A0, eVal, maxIter, tol));

            // (Ford, (2014), Section 18.10)
            // Check NaNs and Infs
            bool hasNan = false;
            label p = 0;
            while (p < eVec.m() && !hasNan)
            {
                const scalar& ev = eVec(p, 0).real();
                if(std::isnan(ev) || std::isinf(ev))
                {
                    hasNan = true;
                }
                ++p;
            }

            // In case of NaNs or Infs, perturb EVal and redo EVec computation
            if (hasNan)
            {
                Random rndGen(1234);
                const scalar perturb
                (
                    (2 + mag(rndGen.GaussNormal<scalar>()))*SMALL
                );

                eVec = eigenvector(A0, eVal, maxIter, tol, perturb);
            }

            //EVecs_.subColumn(k) = Q*eVec;
            for (label i = 0; i < EVecs_.m(); ++i)
            {
                for (label q = 0; q < eVec.m(); ++q)
                {
                    EVecs_(i, k) += Q(i, q)*eVec(q, 0);
                }
            }
        }
    }
    EVecs_.round();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::EigenMatrix<MatrixType>::EigenMatrix()
:
    HessenbergMatrix<MatrixType>
    (
        HessenbergMatrix<MatrixType>::storeMethod::IN_PLACE,
        HessenbergMatrix<MatrixType>::computationQ::FALSE
    ),
    store_(storeMethod::IN_PLACE),
    computeEVecs_(computationEVecs::FALSE),
    H_(HessenbergMatrix<MatrixType>::H_)
{}


template<class MatrixType>
Foam::EigenMatrix<MatrixType>::EigenMatrix
(
    const storeMethod store,
    const computationEVecs computeEVecs
)
:
    HessenbergMatrix<MatrixType>
    (
        static_cast<const typename HessenbergMatrix<MatrixType>::storeMethod>
        (
            store
        ),
        static_cast<const typename HessenbergMatrix<MatrixType>::computationQ>
        (
            computeEVecs
        )
    ),
    store_(store),
    computeEVecs_(computeEVecs),
    H_(HessenbergMatrix<MatrixType>::H_),
    EVals_(2, List<scalar>()),
    EVecs_()
{}


template<class MatrixType>
Foam::EigenMatrix<MatrixType>::EigenMatrix
(
    MatrixType& A,
    const storeMethod store,
    const computationEVecs computeEVecs,
    const label maxIterEVal,
    const scalar tolEVal,
    const label maxIterEVec,
    const scalar tolEVec
)
:
    HessenbergMatrix<MatrixType>
    (
        static_cast<const typename HessenbergMatrix<MatrixType>::storeMethod>
        (
            store
        ),
        static_cast<const typename HessenbergMatrix<MatrixType>::computationQ>
        (
            computeEVecs
        )
    ),
    store_(store),
    computeEVecs_(computeEVecs),
    H_(HessenbergMatrix<MatrixType>::H_),
    EVals_(2, List<scalar>()),
    EVecs_()
{
    decompose(A, maxIterEVal, tolEVal, maxIterEVec, tolEVec);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class MatrixType>
void Foam::EigenMatrix<MatrixType>::decompose
(
    MatrixType& A,
    label maxIterEVal,
    scalar tolEVal,
    label maxIterEVec,
    scalar tolEVec
)
{
    const bool symmetric = A.symmetric();

    if (maxIterEVal < 1)
    {
        // (Kressner, (2005), p. 34)
        maxIterEVal = 30*A.n();
    }

    if (tolEVal < SMALL)
    {
        tolEVal = 1e-5*A.norm();
    }

    // Store a copy of the Hessenberg matrix for eigenvector computations
    MatrixType copyHess;

    // Eigenvalue computations
    HessenbergMatrix<MatrixType>::decompose(A);

    switch (store_)
    {
        case storeMethod::IN_PLACE:
            eigenvalues(A, copyHess, symmetric, maxIterEVal, tolEVal);
            break;

        case storeMethod::OUT_OF_PLACE:
            eigenvalues(H_, copyHess, symmetric, maxIterEVal, tolEVal);
            break;
    }

    // Eigenvector computations
    if (computeEVecs_)
    {
        if (maxIterEVec < 1)
        {
            maxIterEVec = 10;
        }

        if (tolEVec < SMALL)
        {
            tolEVec = 1e-8;
        }

        eigenvectors(copyHess, symmetric, maxIterEVec, tolEVec);
    }
}


template<class MatrixType>
Foam::RectangularMatrix<Foam::scalar> Foam::EigenMatrix<MatrixType>::eigenvector
(
    MatrixType& A0,
    const scalar eVal,
    const label maxIter,
    const scalar tol,
    const scalar perturb
)
{
    Random rndGen(12345);
    const label mRows = A0.m();

    RMatrix change(mRows, 1);
    RMatrix eVec(mRows, 1);
    auto rnd = [&]{ return rndGen.GaussNormal<scalar>(); };
    std::generate
    (
        eVec.begin(),
        eVec.end(),
        rnd
    );
    eVec /= eVec.columnNorm(0);

    //
    for (label i = 0; i < mRows; ++i)
    {
        A0(i, i) -= (eVal + perturb);
    }

    QRMatrix<SMatrix> QRH
    (
        A0,
        QRMatrix<SMatrix>::outputTypes::FULL_QR,
        QRMatrix<SMatrix>::storeMethods::OUT_OF_PLACE,
        QRMatrix<SMatrix>::colPivoting::FALSE
    );
    const SMatrix& R = QRH.R();
    const SMatrix& Q = QRH.Q();

    scalar err = GREAT;
    label i = 0;
    while (tol < err && i < maxIter)
    {
        eVec = Q & eVec;
        eVec = QRH.backSubstitution(R, eVec);
        eVec /= eVec.columnNorm(0);

        //(A0*eVec - eVal*eVec) = (A0 - eVal)*eVec
        change = A0*eVec;
        err = change.columnNorm(0);
        ++i;
    }

    for (label i = 0; i < mRows; ++i)
    {
        A0(i, i) += (eVal + perturb);
    }

    #ifdef FULLDEBUG
    if (i == maxIter)
    {
        Info<< "Not converged for " << i << "th eigenvector" << nl;
    }
    #endif

    return eVec;
}


template<class MatrixType>
Foam::RectangularMatrix<Foam::complex> Foam::EigenMatrix<MatrixType>::
eigenvector
(
    const MatrixType& A0,
    const complex& eVal,
    const label maxIter,
    const scalar tol,
    const scalar perturb
)
{
    Random rndGen(12345);
    const label mRows = A0.m();

    SCMatrix Ac(A0.m());
    auto convertToComplex = [&](const scalar& val) { return complex(val); };
    std::transform
    (
        A0.cbegin(),
        A0.cend(),
        Ac.begin(),
        convertToComplex
    );

    RCMatrix change(mRows, 1);
    RCMatrix eVec(mRows, 1);
    auto rnd = [&]{ return complex(mag(rndGen.GaussNormal<scalar>()), 0); };
    std::generate
    (
        eVec.begin(),
        eVec.end(),
        rnd
    );
    eVec /= complex(eVec.columnNorm(0));

    for (label i = 0; i < mRows; ++i)
    {
        Ac(i, i) -= (eVal + perturb);
    }

    QRMatrix<SCMatrix> QRH
    (
        Ac,
        QRMatrix<SCMatrix>::outputTypes::FULL_QR,
        QRMatrix<SCMatrix>::storeMethods::OUT_OF_PLACE,
        QRMatrix<SCMatrix>::colPivoting::FALSE
    );
    const SCMatrix& R(QRH.R());
    const SCMatrix& Q(QRH.Q());

    scalar err = GREAT;
    label i = 0;
    while (tol < err && i < maxIter)
    {
        eVec = Q & eVec;
        eVec = QRH.backSubstitution(R, eVec);
        eVec /= complex(eVec.columnNorm(0));

        //(A0*eVec - eVal*eVec) = (A0 - eVal)*eVec
        for (label j = 0; j < A0.m(); ++j)
        {
            for (label k = 0; k < eVec.m(); ++k)
            {
                change(j, 0) += A0(j, k)*eVec(k, 0);
            }
        }
        err = change.columnNorm(0);
        ++i;
    }

    #ifdef FULLDEBUG
    if (i == maxIter)
    {
        Info<< "Not converged for " << i << "th eigenvector" << nl;
    }
    #endif

    return eVec;
}


// ************************************************************************* //
