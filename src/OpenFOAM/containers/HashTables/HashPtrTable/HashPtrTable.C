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

#include "error.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable
(
    const HashPtrTable<T, Key, Hash>& ht
)
:
    parent_type(ht.capacity())
{
    for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        const T* ptr = iter.val();
        if (ptr)
        {
            this->set(iter.key(), new T(*ptr));
        }
        else
        {
            this->set(iter.key(), nullptr);
        }
    }
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable
(
    HashPtrTable<T, Key, Hash>&& ht
)
:
    parent_type(std::move(ht))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::~HashPtrTable()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::remove(iterator& iter)
{
    if (iter.good())
    {
        autoPtr<T> aptr(iter.val());
        this->parent_type::erase(iter);
        return aptr;
    }

    return nullptr;
}


template<class T, class Key, class Hash>
Foam::autoPtr<T> Foam::HashPtrTable<T, Key, Hash>::remove(const Key& key)
{
    auto iter = this->find(key);
    return this->remove(iter);
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(iterator& iter)
{
    if (iter.good())
    {
        T* ptr = iter.val();

        if (this->parent_type::erase(iter))
        {
            if (ptr)
            {
                delete ptr;
            }

            return true;
        }
    }

    return false;
}


template<class T, class Key, class Hash>
bool Foam::HashPtrTable<T, Key, Hash>::erase(const Key& key)
{
    auto iter = this->find(key);
    return this->erase(iter);
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::clear()
{
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        delete iter.val();
    }

    this->parent_type::clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::operator=
(
    const HashPtrTable<T, Key, Hash>& rhs
)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted copy assignment to self"
            << abort(FatalError);
    }

    this->clear();

    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        const T* ptr = iter.val();
        if (ptr)
        {
            this->set(iter.key(), new T(*ptr));
        }
        else
        {
            this->set(iter.key(), nullptr);
        }
    }
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::operator=
(
    HashPtrTable<T, Key, Hash>&& rhs
)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted move assignment to self"
            << abort(FatalError);
    }

    this->clear();
    this->transfer(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashPtrTableIO.C"

// ************************************************************************* //
