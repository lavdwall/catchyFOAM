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

inline Foam::pdepArrheniusReactionRate::pdepArrheniusReactionRate
(
    const scalarList p,
    const scalarList A,
    const scalarList beta,
    const scalarList Ta
)
:
    p_(p),
    A_(A),
    beta_(beta),
    Ta_(Ta)
{}


inline Foam::pdepArrheniusReactionRate::pdepArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    p_(dict.lookup("p")),
    A_(dict.lookup("A")),
    beta_(dict.lookup("beta")),
    Ta_(dict.lookup("Ta"))
{
    if( (p_.size() != A_.size()) || (p_.size() != beta_.size()) || (p_.size() != Ta_.size()) )
    {        
        FatalErrorInFunction
            << "Error in constructing pdepArrheniusReactionRate object: arrays not the same size!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pdepArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    label mini = 0; 
    scalar mindeltap = 1e15;
    forAll(p_, i)
    {
        if (mag(p - p_[i]) < mag(mindeltap))
        {
            mini = i;
            mindeltap = p - p_[i];
        }      
    }

    scalar aki = A_[mini];
    if (mag(beta_[mini]) > vSmall)
    {
        aki *= pow(T, beta_[mini]);
    }
    if (mag(Ta_[mini]) > vSmall)
    {
        aki *= exp(-Ta_[mini]/T);
    }

    if ( mini==0 || mini==p_.size()-1)
    {
        return aki;
    }
    else
    {
        label j = (mindeltap>0 ? mini+1 : mini-1);

        scalar akj = A_[j];
        if (mag(beta_[j]) > vSmall)
        {
            akj *= pow(T, beta_[j]);
        }
        if (mag(Ta_[j]) > vSmall)
        {
            akj *= exp(-Ta_[j]/T);
        }
        
        if (mindeltap>0) {
            scalar ln_ak = log(aki) + (log(akj)-log(aki))*(log(p)-log(p_[mini]))/(log(p_[j])-log(p_[mini]));
            return exp(ln_ak);
        } else {
            scalar ln_ak = log(akj) + (log(aki)-log(akj))*(log(p)-log(p_[j]))/(log(p_[mini])-log(p_[j]));
            return exp(ln_ak);
        }
    }
}


inline Foam::scalar Foam::pdepArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    label mini = 0; 
    scalar mindeltap = 1e15;
    forAll(p_, i)
    {
        if (mag(p - p_[i]) < mag(mindeltap))
        {
            mini = i;
            mindeltap = p - p_[i];
        }        
    }

    scalar aki = A_[mini];
    if (mag(beta_[mini]) > vSmall)
    {
        aki *= pow(T, beta_[mini]);
    }
    if (mag(Ta_[mini]) > vSmall)
    {
        aki *= exp(-Ta_[mini]/T);
    }
    scalar daki = aki*(beta_[mini]+Ta_[mini]/T)/T;

    if ( mini==0 || mini==p_.size()-1)
    {
        return daki;
    }
    else
    {
        label j = (mindeltap>0 ? mini+1 : mini-1);

        scalar akj = A_[j];
        if (mag(beta_[j]) > vSmall)
        {
            akj *= pow(T, beta_[j]);
        }
        if (mag(Ta_[j]) > vSmall)
        {
            akj *= exp(-Ta_[j]/T);
        }
        scalar dakj = akj*(beta_[j]+Ta_[j]/T)/T;
        
        scalar ln_ak = log(aki) + (log(akj)-log(aki))*(log(p)-log(p_[mini]))/(log(p_[j])-log(p_[mini]));
        return exp(ln_ak)*(1.0/aki*daki + (log(p)-log(p_[mini]))/(log(p_[j])-log(p_[mini]))*(1.0/akj*dakj - 1.0/aki*daki));
    }
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::pdepArrheniusReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::pdepArrheniusReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::pdepArrheniusReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label
) const
{
    return 0;
}


inline void Foam::pdepArrheniusReactionRate::write(Ostream& os) const
{
    os.writeKeyword("p") << p_ << token::END_STATEMENT << nl;
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pdepArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
