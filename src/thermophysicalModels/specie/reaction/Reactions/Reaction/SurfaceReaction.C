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

#include "SurfaceReaction.H"
#include "DynamicList.H"
#include "IOobject.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::SurfaceReaction<ReactionThermo>::SurfaceReaction
(
    const Reaction<ReactionThermo>& reaction,
    const speciesTable& gases,
    const List<specieCoeffs>& glhs,
    const List<specieCoeffs>& grhs
)
:
    Reaction<ReactionThermo>(reaction),
    gases_(gases),
    nGases_(gases_.size()),
    glhs_(glhs),
    grhs_(grhs)
{}


template<class ReactionThermo>
Foam::SurfaceReaction<ReactionThermo>::SurfaceReaction
(
    const SurfaceReaction<ReactionThermo>& r,
    const speciesTable& gases
)
:
    Reaction<ReactionThermo>(r),
    gases_(gases),
    nGases_(gases_.size()),
    glhs_(r.glhs_),
    grhs_(r.grhs_)
{}


template<class ReactionThermo>
Foam::SurfaceReaction<ReactionThermo>::SurfaceReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>
    (
        dict.dictName(),
        species,
        dict.lookupOrDefault<scalar>("Tlow", Reaction<ReactionThermo>::TlowDefault),
        dict.lookupOrDefault<scalar>("Thigh", Reaction<ReactionThermo>::ThighDefault),
        thermoDatabase,
        dict,
        false
    ),
    gases_(dictionary(IFstream(fileName(dict.parent().parent().lookup("gasesThermoLocation")).expand())()).subDict("gas").lookup("species")),
    nGases_(gases_.size())
{
    List<specieCoeffs> lhstemp_;
    List<specieCoeffs> rhstemp_;
    specieCoeffs::setLRhsGS
    (
        IStringStream(dict.lookup("reaction"))(),
        gases_,
        this->species(),
        lhstemp_,
        rhstemp_
    );
    forAll(lhstemp_,i)
    {
        if(lhstemp_[i].isgas == true) glhs_.append(lhstemp_[i]);
        else this->lhs_.append(lhstemp_[i]);
    }
    forAll(rhstemp_,i)
    {
        if(rhstemp_[i].isgas == true) grhs_.append(rhstemp_[i]);
        else this->rhs_.append(rhstemp_[i]);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::autoPtr<Foam::SurfaceReaction<ReactionThermo>>
Foam::SurfaceReaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
{
    const word& reactionTypeName = dict.lookup("type");

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. Reaction names used to have "Reaction"
    // (Reaction<ReactionThermo>::typeName_()) appended. This was removed as it
    // is unnecessary given the context in which the reaction is specified. If
    // this reaction name was not found, search also for the old name.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<SurfaceReaction<ReactionThermo>>
    (
        cstrIter()(species, thermoDatabase, dict)
    );
}

template<class ReactionThermo>
Foam::autoPtr<Foam::SurfaceReaction<ReactionThermo>>
Foam::SurfaceReaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // If the objectRegistry constructor table is empty
    // use the dictionary constructor table only
    if (!objectRegistryConstructorTablePtr_)
    {
        return New(species, thermoDatabase, dict);
    }

    const word& reactionTypeName = dict.lookup("type");

    typename objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. See above.
    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        cstrIter = objectRegistryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(reactionTypeName);

        // Backwards compatibility check. See above.
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            cstrIter = dictionaryConstructorTablePtr_->find
            (
                reactionTypeName.removeTrailing(typeName_())
            );
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown reaction type "
                << reactionTypeName << nl << nl
                << "Valid reaction types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << objectRegistryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<SurfaceReaction<ReactionThermo>>
        (
            cstrIter()(species, thermoDatabase, dict)
        );
    }

    return autoPtr<SurfaceReaction<ReactionThermo>>
    (
        cstrIter()(species, thermoDatabase, ob, dict)
    );
}

template<class ReactionThermo>
Foam::autoPtr<Foam::SurfaceReaction<ReactionThermo>>
Foam::SurfaceReaction<ReactionThermo>::New
(
    const speciesTable& species,
    const PtrList<ReactionThermo>& speciesThermo,
    const dictionary& dict
)
{
    HashPtrTable<ReactionThermo> thermoDatabase;
    forAll(speciesThermo, i)
    {
        thermoDatabase.insert
        (
            speciesThermo[i].name(),
            speciesThermo[i].clone().ptr()
        );
    }

    return New(species, thermoDatabase, dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
const Foam::List<Foam::specieCoeffs>&
Foam::SurfaceReaction<ReactionThermo>::glhs() const
{
    return glhs_;
}


template<class ReactionThermo>
const Foam::List<Foam::specieCoeffs>&
Foam::SurfaceReaction<ReactionThermo>::grhs() const
{
    return grhs_;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::SurfaceReaction<ReactionThermo>::
gasSpecies() const
{
    return gases_;
}

template<class ReactionThermo>
void Foam::SurfaceReaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalar omegai = omega
    (
        p, T, c, li, pf, cf, lRef, pr, cr, rRef
    );

    forAll(this->lhs(), s)
    {
        label si = this->lhs()[s].index + nGases_;
        scalar sl = this->lhs()[s].stoichCoeff;
        dcdt[si] -= sl*omegai; // [kmol/m2]
    }
    forAll(this->glhs(), g)
    {
        label gi = this->glhs()[g].index;
        scalar gl = this->glhs()[g].stoichCoeff;
        dcdt[gi] -= gl*omegai; // [kmol/m2]
    }
    forAll(this->rhs(), s)
    {
        label si = this->rhs()[s].index + nGases_;
        scalar sr = this->rhs()[s].stoichCoeff;
        dcdt[si] += sr*omegai; // [kmol/m2]
    }
    forAll(this->grhs(), g)
    {
        label gi = this->grhs()[g].index;
        scalar gr = this->grhs()[g].stoichCoeff;
        dcdt[gi] += gr*omegai; // [kmol/m2]
    }
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalar clippedT = min(max(T, this->Tlow()), this->Thigh());

    const scalar kf = this->kf(p, clippedT, c, li);
    const scalar kr = this->kr(kf, p, clippedT, c, li);

    pf = 1;
    pr = 1;

    const label Nl = this->lhs().size();
    const label Nr = this->rhs().size();
    const label Ngl = this->glhs().size();
    const label Ngr = this->grhs().size();

    label slRef = 0;
    lRef = this->lhs()[slRef].index + nGases_;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = this->lhs()[s].index + nGases_;

        if (c[si] < c[lRef])
        {
            const scalar exp = this->lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = this->lhs()[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    for (label g = 0; g < Ngl; g++)
    {
        label gi = this->glhs()[g].index;
        const scalar exp = this->glhs()[g].exponent;
        pf *= pow(max(c[gi], 0), exp);
    }
    cf = max(c[lRef], 0);
    {
        const scalar exp = this->lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = this->rhs()[srRef].index + nGases_;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = this->rhs()[s].index + nGases_;
        if (c[si] < c[rRef])
        {
            const scalar exp = this->rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = this->rhs()[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    for (label g = 0; g < Ngr; g++)
    {
        label gi = this->grhs()[g].index;
        const scalar exp = this->grhs()[g].exponent;
        pr *= pow(max(c[gi], 0), exp);
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = this->rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo>
void Foam::SurfaceReaction<ReactionThermo>::dwdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarSquareMatrix& J,
    scalarField& dcdt,
    scalar& omegaI,
    scalar& kfwd,
    scalar& kbwd,
    const bool reduced,
    const List<label>& c2s
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    omegaI = omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef);

    forAll(this->lhs(), s)
    {
        label si = this->lhs()[s].index + nGases_;
        scalar sl = this->lhs()[s].stoichCoeff;
        dcdt[si] -= sl*omegaI; // [kmol/m2]
    }
    forAll(this->glhs(), g)
    {
        label gi = this->glhs()[g].index;
        scalar gl = this->glhs()[g].stoichCoeff;
        dcdt[gi] -= gl*omegaI; // [kmol/m2]
    }
    forAll(this->rhs(), s)
    {
        label si = this->rhs()[s].index + nGases_;
        scalar sr = this->rhs()[s].stoichCoeff;
        dcdt[si] += sr*omegaI; // [kmol/m2]
    }
    forAll(this->grhs(), g)
    {
        label gi = this->grhs()[g].index;
        scalar gr = this->grhs()[g].stoichCoeff;
        dcdt[gi] += gr*omegaI; // [kmol/m2]
    }

    kfwd = this->kf(p, T, c, li);
    kbwd = this->kr(kfwd, p, T, c, li);

    forAll(this->lhs(), j)
    {
        const label sj = this->lhs()[j].index + nGases_;
        scalar kf = kfwd;
        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar el = this->lhs()[i].exponent;
            if (i == j)
            {
                if (el < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kf *= el*pow(c[si] + VSMALL, el - 1);
                    }
                    else
                    {
                        kf = 0;
                    }
                }
                else
                {
                    kf *= el*pow(c[si], el - 1);
                }
            }
            else
            {
                kf *= pow(c[si], el);
            }
        }
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar el = this->glhs()[i].exponent;
            kf *= pow(c[si], el);
        }

        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar sl = this->lhs()[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar sl = this->glhs()[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar sr = this->rhs()[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar sr = this->grhs()[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }

    forAll(this->glhs(), j)
    {
        const label sj = this->glhs()[j].index;
        scalar kf = kfwd;
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar el = this->glhs()[i].exponent;
            if (i == j)
            {
                if (el < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kf *= el*pow(c[si] + VSMALL, el - 1);
                    }
                    else
                    {
                        kf = 0;
                    }
                }
                else
                {
                    kf *= el*pow(c[si], el - 1);
                }
            }
            else
            {
                kf *= pow(c[si], el);
            }
        }
        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar el = this->lhs()[i].exponent;
            kf *= pow(c[si], el);
        }

        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar sl = this->lhs()[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar sl = this->glhs()[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar sr = this->rhs()[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar sr = this->grhs()[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }

    forAll(this->rhs(), j)
    {
        const label sj = this->rhs()[j].index + nGases_;
        scalar kr = kbwd;
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar er = this->rhs()[i].exponent;
            if (i == j)
            {
                if (er < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kr *= er*pow(c[si] + VSMALL, er - 1);
                    }
                    else
                    {
                        kr = 0;
                    }
                }
                else
                {
                    kr *= er*pow(c[si], er - 1);
                }
            }
            else
            {
                kr *= pow(c[si], er);
            }
        }
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar er = this->grhs()[i].exponent;
            kr *= pow(c[si], er);
        }

        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar sl = this->lhs()[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar sl = this->glhs()[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar sr = this->rhs()[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar sr = this->grhs()[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
    }

    forAll(this->grhs(), j)
    {
        const label sj = this->grhs()[j].index;
        scalar kr = kbwd;
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar er = this->grhs()[i].exponent;
            if (i == j)
            {
                if (er < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kr *= er*pow(c[si] + VSMALL, er - 1);
                    }
                    else
                    {
                        kr = 0;
                    }
                }
                else
                {
                    kr *= er*pow(c[si], er - 1);
                }
            }
            else
            {
                kr *= pow(c[si], er);
            }
        }
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar er = this->rhs()[i].exponent;
            kr *= pow(c[si], er);
        }

        forAll(this->lhs(), i)
        {
            const label si = this->lhs()[i].index + nGases_;
            const scalar sl = this->lhs()[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(this->glhs(), i)
        {
            const label si = this->glhs()[i].index;
            const scalar sl = this->glhs()[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(this->rhs(), i)
        {
            const label si = this->rhs()[i].index + nGases_;
            const scalar sr = this->rhs()[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
        forAll(this->grhs(), i)
        {
            const label si = this->grhs()[i].index;
            const scalar sr = this->grhs()[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
    }

    // When coverage dependencies are involved, additional terms are added
    // beta function returns an empty list when coverage dependencies are not involved
    const List<Tuple2<label, scalarField>>& covdep = this->covdep();
    if (notNull(covdep))
    {
        // This temporary array needs to be cached for efficiency
        scalarField dcidc(c.size(), 0.0);
        this->dcidc(p, T, c, li, dcidc);

        forAll(covdep, j)
        {
            label sj = covdep[j].first();
            if (sj != -1)
            {
                forAll(this->lhs(), i)
                {
                    const label si = this->lhs()[i].index + nGases_;
                    const scalar sl = this->lhs()[i].stoichCoeff;
                    J(si, sj) -= sl*dcidc[sj]*omegaI;
                }
                forAll(this->glhs(), i)
                {
                    const label si = this->glhs()[i].index;
                    const scalar sl = this->glhs()[i].stoichCoeff;
                    J(si, sj) -= sl*dcidc[sj]*omegaI;
                }
                forAll(this->rhs(), i)
                {
                    const label si = this->rhs()[i].index + nGases_;
                    const scalar sr = this->rhs()[i].stoichCoeff;
                    J(si, sj) += sr*dcidc[sj]*omegaI;
                }
                forAll(this->grhs(), i)
                {
                    const label si = this->grhs()[i].index;
                    const scalar sr = this->grhs()[i].stoichCoeff;
                    J(si, sj) += sr*dcidc[sj]*omegaI;
                }
            }
        }
    }
}


template<class ReactionThermo>
void Foam::SurfaceReaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    writeEntry
    (
        os,
        "reaction",
        solidReactionStr(reaction)
    );
}


template<class ReactionThermo>
Foam::string Foam::SurfaceReaction<ReactionThermo>::solidReactionStr
(
    OStringStream& reaction
) const
{
    specieCoeffs::reactionStr(reaction, this->species(), this->lhs());
    if (glhs().size() > 0)
    {
        reaction << " + ";
        specieCoeffs::reactionStr(reaction, this->gasSpecies(), this->glhs());
    }
    reaction << " = ";
    specieCoeffs::reactionStr(reaction, this->species(), this->rhs());
    if (grhs().size() > 0)
    {
        reaction << " + ";
        specieCoeffs::reactionStr(reaction, this->gasSpecies(), this->grhs());
    }
    return reaction.str();

}

// ************************************************************************* //
