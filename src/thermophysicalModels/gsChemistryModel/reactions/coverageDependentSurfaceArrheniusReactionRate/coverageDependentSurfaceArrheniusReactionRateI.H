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

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::coverageDependentSurfaceArrheniusReactionRate::coverageDependentSurfaceArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const speciesTable& species,
    const scalar sden,
    const label nGases,
    const List<Tuple2<word, scalarField>> cov
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    species_(species),
    sden_(sden),
    nGases_(nGases)
{
    forAll(cov, i)
    {
        cov_.append
        (
            Tuple2<label, scalarField>
            (
                species[cov[i].first()]+nGases_,
                cov[i].second()
            )
        );
    }
}


inline Foam::coverageDependentSurfaceArrheniusReactionRate::coverageDependentSurfaceArrheniusReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ta_(readScalar(dict.lookup("Ta"))),
    species_(species),
    sden_(readScalar(dict.lookup("sden"))),
    nGases_(readLabel(dict.lookup("ngasOffset")))
{
    if (dict.found("coverages"))
    {
        List<Tuple2<word, scalarField>> coeffs(dict.lookup("coverages"));

        forAll(coeffs, i)
        {
            cov_.append
            (
                Tuple2<label, scalarField>
                (
                    species_[coeffs[i].first()]+nGases_,
                    coeffs[i].second()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::coverageDependentSurfaceArrheniusReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    forAll(cov_, i)
    {
        scalar covi = c[cov_[i].first()]/sden_*cov_[i].second()[3];
        scalar ai = cov_[i].second()[0];
        scalar mi = cov_[i].second()[1];
        scalar Ti = cov_[i].second()[2];

        ak *= pow(10.0, ai*covi)*pow(covi, mi)*exp(-Ti*covi/T);
    }

    return ak;
}


inline Foam::scalar Foam::coverageDependentSurfaceArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    scalar covprodk = 1.0;
    scalar covsumk = 0.0;
    forAll(cov_, i)
    {
        scalar covi = c[cov_[i].first()]/sden_*cov_[i].second()[3];
        scalar ai = cov_[i].second()[0];
        scalar mi = cov_[i].second()[1];
        scalar Ti = cov_[i].second()[2];

        covprodk *= pow(10.0, ai*covi)*pow(covi, mi)*exp(-Ti*covi/T);
        covsumk += Ti*covi/pow(T, 2.0);
    }

    return ak*(beta_+Ta_/T)/T*covprodk + ak*covprodk*covsumk;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::coverageDependentSurfaceArrheniusReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::coverageDependentSurfaceArrheniusReactionRate::covdep() const
{
    return cov_;
}


inline void Foam::coverageDependentSurfaceArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    forAll(cov_, i)
    {
        scalar covi = c[cov_[i].first()]/sden_*cov_[i].second()[3];
        scalar ai = cov_[i].second()[0];
        scalar mi = cov_[i].second()[1];
        scalar Ti = cov_[i].second()[2];

        if (mag(covi) > vSmall)
            dcidc[cov_[i].first()] = (log(10.0)*T*ai*covi+T*mi-Ti*covi)/(covi*T)*cov_[i].second()[3]/sden_;
    }
}


inline Foam::scalar Foam::coverageDependentSurfaceArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    scalar covsumk = 0.0;
    forAll(cov_, i)
    {
        scalar covi = c[cov_[i].first()]/sden_*cov_[i].second()[3];
        scalar Ti = cov_[i].second()[2];

        covsumk += Ti*covi/pow(T, 2.0);
    }
    return covsumk;
}


inline void Foam::coverageDependentSurfaceArrheniusReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
    os.writeKeyword("sden") << sden_ << token::END_STATEMENT << nl;
    os.writeKeyword("ngasOffset") << nGases_ << token::END_STATEMENT << nl;

    List<Tuple2<word, scalarField>> coeffs(cov_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[cov_[i].first()-nGases_];
        coeffs[i].second() = cov_[i].second();
    }
    os.writeKeyword("coverages") << coeffs << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const coverageDependentSurfaceArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
