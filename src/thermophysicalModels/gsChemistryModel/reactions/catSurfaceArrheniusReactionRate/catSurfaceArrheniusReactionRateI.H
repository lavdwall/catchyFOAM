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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::catSurfaceArrheniusReactionRate::catSurfaceArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta)
{}


inline Foam::catSurfaceArrheniusReactionRate::catSurfaceArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(readScalar(dict.lookup("A"))),
    beta_(readScalar(dict.lookup("beta"))),
    Ta_(readScalar(dict.lookup("Ta")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::catSurfaceArrheniusReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField&,
    const label
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

    return ak;
}


inline Foam::scalar Foam::catSurfaceArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
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

    return ak*(beta_+Ta_/T)/T;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::catSurfaceArrheniusReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::catSurfaceArrheniusReactionRate::covdep() const
{
    return NullObjectRef<List<Tuple2<label, scalarField>>>();
}


inline void Foam::catSurfaceArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::catSurfaceArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label
) const
{
    return 0;
}


inline void Foam::catSurfaceArrheniusReactionRate::write(Ostream& os) const
{
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const catSurfaceArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
