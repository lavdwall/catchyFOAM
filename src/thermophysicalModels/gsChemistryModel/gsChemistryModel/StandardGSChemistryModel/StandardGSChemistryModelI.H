/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::label
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RRg()
{
    return RRg_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RRs()
{
    return RRs_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<Foam::Reaction<GThermoType>>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::gasReactions() const
{
    return gasReactions_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<Foam::SurfaceReaction<SThermoType>>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solidReactions() const
{
    return solidReactions_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<GThermoType>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::gasSpecieThermo() const
{
    return gasSpecieThermo_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<SThermoType>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solidSpecieThermo() const
{
    return solidSpecieThermo_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::label
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::label
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::nReaction() const
{
    return nReaction_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalar
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Treact() const
{
    return Treact_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalar&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Treact()
{
    return Treact_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RRg
(
    const label i
) const
{
    return RRg_[i];
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RRs
(
    const label i
) const
{
    return RRs_[i];
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RR
(
    const label i
) const
{
    if (i<nGasSpecie_)
    {
        return RRg(i);
    }
    else
    {
        return RRs(i-nGasSpecie_);
    }
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::RR
(
    const label i
)
{
    if (i<nGasSpecie_)
    {
        return RRg_[i];
    }
    else
    {
        return RRs_[i];
    }
}


// ************************************************************************* //
