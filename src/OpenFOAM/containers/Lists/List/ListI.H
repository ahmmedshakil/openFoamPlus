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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::doAlloc()
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];
    }
}


template<class T>
inline void Foam::List<T>::reAlloc(const label len)
{
    if (this->size_ != len)
    {
        clear();
        this->size_ = len;
        doAlloc();
    }
}


template<class T>
template<class List2>
inline void Foam::List<T>::copyList(const List2& list)
{
    const label len = this->size_;

    for (label i=0; i<len; ++i)
    {
        this->operator[](i) = list[i];
    }
}


template<class T>
template<class InputIterator>
inline Foam::List<T>::List
(
    InputIterator begIter,
    InputIterator endIter,
    const label len
)
:
    UList<T>(nullptr, len)
{
    if (this->size_)
    {
        doAlloc();

        InputIterator iter = begIter;
        for (label i = 0; i < len; ++i)
        {
            this->operator[](i) = *iter;
            ++iter;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::List<T>::List() noexcept
{}


template<class T>
inline Foam::autoPtr<Foam::List<T>> Foam::List<T>::clone() const
{
    return autoPtr<List<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::List<T>& Foam::List<T>::null()
{
    return NullObjectRef<List<T>>();
}


template<class T>
inline void Foam::List<T>::clear()
{
    if (this->v_)
    {
        delete[] this->v_;
        this->v_ = nullptr;
    }

    this->size_ = 0;
}


namespace Foam
{
    // Template specialization for bool. Fills with false
    template<>
    inline void List<bool>::resize(const label newSize)
    {
        this->resize(newSize, false);
    }
}


template<class T>
inline void Foam::List<T>::resize(const label len)
{
    this->doResize(len);
}


template<class T>
inline void Foam::List<T>::setSize(const label len)
{
    this->resize(len);
}


template<class T>
inline void Foam::List<T>::setSize(const label len, const T& val)
{
    this->resize(len, val);
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    const label n = this->size();

    if (i >= n)
    {
        resize(2*n);
    }

    return UList<T>::operator[](i);
}


template<class T>
inline void Foam::List<T>::append(const T& val)
{
    resize(this->size() + 1, val);  // copy element
}


template<class T>
inline void Foam::List<T>::append(T&& val)
{
    const label idx = this->size();
    resize(idx + 1);

    this->operator[](idx) = std::move(val);  // move assign element
}


template<class T>
inline void Foam::List<T>::append(const UList<T>& list)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "attempted appending to self" << abort(FatalError);
    }

    label idx = this->size();
    const label n = list.size();

    resize(idx + n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = list[i];  // copy element
    }
}


template<class T>
template<class Addr>
inline void Foam::List<T>::append(const IndirectListBase<T, Addr>& list)
{
    label idx = this->size();
    const label n = list.size();

    resize(idx + n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = list[i];  // copy element
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::operator=(const zero)
{
    UList<T>::operator=(Zero);
}


// ************************************************************************* //
