/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "ParticleGSChemistryModel.H"
#include "multiComponentGSMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"


// * * * * * * * * * * *  Private member functions  * * * * * * * * * * * * * //
template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
initializeSurface()
{
    initialize_ = true;
    if (this->lookupOrDefault("initializeCoveragesQSSA", false))
    {
        Info << "Initializing coverages using QSSA" << endl;
        const dictionary& QSSADict = this->subDict("coveragesQSSA");
        scalar treshold_ = QSSADict.lookupOrDefault("treshold", 1e-3);
        scalar deltaT_ = QSSADict.lookupOrDefault("deltaT", this->mesh().time().deltaTValue());
        label maxIter_ = QSSADict.lookupOrDefault("maxIter", 1000);

        if (QSSADict.lookupOrDefault("uniformCoverages", false))
        {
            Info << "Setting uniform coverages" << endl;
            label cellC = 0;
            solveSurfaceQSSA(cellC, treshold_, maxIter_, deltaT_);

            forAll(Ys_, i)
            {
                forAll(Ys_[i], celli)
                {
                    Ys_[i][celli] = Ys_[i][cellC];
                }
                forAll(Ys_[i].boundaryField(), patchi)
                {
                    fvPatchScalarField& Ypi = Ys_[i].boundaryFieldRef()[patchi];
                    forAll(Ypi, facei)
                    {
                        Ypi[facei] = Ypi.patchInternalField()()[facei];
                    }
                }
            }
        }
        else
        {
            forAll(Ys_[0], celli)
            {
                solveSurfaceQSSA(celli, treshold_, maxIter_, deltaT_);
            }
            forAll(Ys_, i)
            {
                forAll(Ys_[i].boundaryField(), patchi)
                {
                    fvPatchScalarField& Ypi = Ys_[i].boundaryFieldRef()[patchi];
                    forAll(Ypi, facei)
                    {
                        Ypi[facei] = Ypi.patchInternalField()()[facei];
                    }
                }
            }
        }
    }
    initialize_ = false;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
solveSurfaceQSSA
(
    const label celli,
    const scalar treshold,
    const label maxiter,
    const scalar deltaT
)
{
    if (!this->chemistry_)
    {
        return;
    }

    scalar Ti = this->thermo().T()[celli];
    scalar pi = this->thermo().p()[celli];
    scalar rhoi =
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).gasCellMixture(celli).rho(pi, Ti);

    scalarField c_old(nSpecie_);

    for (label i=0; i<nGasSpecie_; i++)
    {
        const scalar Yi = Y_[i][celli];
        c_[i] = rhoi*Yi/gasSpecieThermo_[i].W();
        c_old[i] = c_[i];
    }
    for (label i=0; i<this->nSolidSpecie_; i++)
    {
        const scalar Ysi = Ys_[i][celli];
        c_[i + nGasSpecie_] = Ysi*solidSpecieThermo_[i].sden()/solidSpecieThermo_[i].size();
        c_old[i + nGasSpecie_] = c_[i + nGasSpecie_];
    }

    label iter = 0;
    scalar deltamax = 1000;
    while ( (deltamax > treshold) && (iter < maxiter))
    {
        // Initialise time progress
        scalar timeLeft = deltaT;

        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
            timeLeft -= dt;
        }

        deltamax = 0.0;
        for (label i=0; i<this->nSolidSpecie_; i++)
        {
            if (c_[i + nGasSpecie_] != 0.0)
            {
                deltamax = max(deltamax, mag((c_[i + nGasSpecie_]-c_old[i + nGasSpecie_])/c_[i + nGasSpecie_]));
            }
            c_old[i + nGasSpecie_] = c_[i + nGasSpecie_];
        }

        iter++;
    }

    for (label i=0; i<this->nSolidSpecie_; i++)
    {
        Ys_[i][celli] = c_[i+nGasSpecie_]/solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::ParticleGSChemistryModel
(
    rhoReactionThermo& thermo
)
:
    ReactionThermo(thermo),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    Ys_(this->thermo().composition().Ys()),
    gasSpecieThermo_
    (
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
            (this->thermo()).specieThermos()
    ),
    solidSpecieThermo_
    (
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
            (this->thermo()).solidSpecieThermos()
    ),
    gasReactions_
    (
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).species(),
        gasSpecieThermo_,
        this->mesh(),
        this->subDict("gas")
    ),
    solidReactions_
    (
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).solidSpecies(),
        solidSpecieThermo_,
        this->mesh(),
        this->subDict("solid")
    ),

    nGasSpecie_(this->thermo().composition().Y().size()),
    nSolidSpecie_(this->thermo().composition().Ys().size()),
    nSpecie_(nGasSpecie_+nSolidSpecie_),
    nGasReaction_(gasReactions_.size()),
    nSolidReaction_(solidReactions_.size()),
    nReaction_(nGasReaction_+nSolidReaction_),
    Treact_
    (
        ReactionThermo::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    c_(nSpecie_+2),
    dcdt_(nSpecie_+2),
    initialize_(false),
    catarea_(this->lookupOrDefault("catalyticSurfaceAreaPerVolume", 1.0)),
    por_(this->lookupOrDefault("catalystPorosity", 0.0)),
    rhogs_(0.0),
    instantaneousRates_(this->lookupOrDefault("instantaneousRates", false))
{
    Info<< "ParticleGSChemistryModel: " << nl
        << indent << "Number of gas species = " << nGasSpecie_ << nl
        << indent << "Number of solid species = " << nSolidSpecie_ << nl
        << indent << "Number of gas reactions = " << nGasReaction_ << nl
        << indent << "Number of solid reactions = " << nSolidReaction_ << nl
        << indent << "Catalyst porosity = " << por_ << nl
        << indent << "Catalyst surface area to volume = " << catarea_ << endl;

    if (instantaneousRates_)
        Info << "Using instantaneous reaction rates" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
~ParticleGSChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{

    dcdt = Zero;

    forAll(solidReactions_, i)
    {
        const SurfaceReaction<SThermoType>& R = solidReactions_[i];

        R.omega(p, T, c, li, dcdt);
    }

    for (label i=0; i<nGasSpecie_; i++)
    {
        dcdt[i] *= catarea_*(initialize_? 0.0 : 1.0);
    }

    if (por_*(initialize_? 0.0 : 1.0) > 0.0)
    {
        scalarField dcdtg(nGasSpecie_, 0.0);

        forAll(gasReactions_, i)
        {
            const Reaction<GThermoType>& R = gasReactions_[i];

            R.omega(p, T, c, li, dcdtg);
        }

        forAll(dcdtg, i)
        {
            dcdt[i] += dcdtg[i]*por_;
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;

    omega(p, T, c_, li, dcdt);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    J = Zero;
    dcdt = Zero;

    scalar omegaI = 0;
    List<label> dummy;
    forAll(solidReactions_, sri)
    {
        const SurfaceReaction<SThermoType>& R = solidReactions_[sri];
        scalar kfwd, kbwd;
        R.dwdc(p, T, c_, li, J, dcdt, omegaI, kfwd, kbwd, false, dummy);
    }
    for (label j=0; j<nGasSpecie_; j++)
    {
        forAll(c_, i)
        {
            J(j,i) *= catarea_;
        }
        dcdt[j] *= catarea_*(initialize_? 0.0 : 1.0);
    }
    if (por_*(initialize_? 0.0 : 1.0)>0.0)
    {
        scalarSquareMatrix Jg(nGasSpecie_+1, 0.0);
        scalarField dcdtg(nGasSpecie_+1, 0.0);

        forAll(gasReactions_, gri)
        {
            const Reaction<GThermoType>& R = gasReactions_[gri];
            scalar kfwd, kbwd;
            R.dwdc(p, T, c_, li, Jg, dcdtg, omegaI, kfwd, kbwd, false, dummy);
        }

        for (label j=0; j<nGasSpecie_; j++)
        {
            for (label i=0; i<nGasSpecie_; i++)
            {
                J(j,i) += Jg(j,i)*por_;
            }
            dcdt[j] += dcdtg[j]*por_;
            J(j,nSpecie_) += Jg(j,nGasSpecie_)*por_;
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::tc() const
{
    NotImplemented;

    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    return ttc;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Qdot() const
{
    NotImplemented;
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );
    return tQdot;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    NotImplemented;

    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    return tRR;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::calculate()
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
template<class DeltaTType>
Foam::scalar Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    NotImplemented;
    scalar deltaTMin = great;
    return deltaTMin;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setGasMassFractions
(
    const PtrList<volScalarField>& Yg
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
disableSurfaceSpeciesEqn()
{
    Info << "Disabling solution of surface species equation" << nl
        << indent << "NOTE: this means that values of surface species fractions "
        << "will be assigned after solution of the chemistry" << endl;
    for (label i=0; i<nSolidSpecie_; i++)
    {
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).setInactive(i+nGasSpecie_);
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
correctRates
(
    PtrList<volScalarField>& Yg,
    const volScalarField& rhog,
    const volScalarField& Tg,
    const volScalarField& pg
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
correctCatalyticWallFluxes
(
    PtrList<volScalarField>& Yg,
    const volScalarField& rhog,
    volScalarField& heg,
    volScalarField& Tg,
    const volScalarField& pg,
    const wordList& patches
)
{
   NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallSpeciesFlux
(
    const word speciei,
    const label patchi
) const
{
    NotImplemented;
    return c_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallHeatFlux
(
    const label patchi
) const
{
    NotImplemented;
    return c_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setAlpha
(
    const volScalarField& alpha
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTemperature
(
    const volScalarField& Tg
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::ParticleGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
getRatesQdotI
(
    const scalar deltaT,
    const label celli,
    scalarField& RR,
    const scalarField& Yg,
    scalarField& Ys,
    const scalar rhog,
    const scalar Tg,
    scalar Ts,
    scalar p,
    const bool ini
)
{
    // Checks
    if (deltaT < 0.0)
      FatalErrorInFunction << "Negative deltaT" << abort(FatalError);

    if (celli < 0 || celli>=this->mesh().V().size())
      FatalErrorInFunction << "Invalid celli" << abort(FatalError);

    if (RR.size() != nGasSpecie_)
      FatalErrorInFunction << "Wrong RR size" << abort(FatalError);

    if (Yg.size() != nGasSpecie_)
      FatalErrorInFunction << "Wrong Yg size" << abort(FatalError);

    if (Ys.size() != nSolidSpecie_)
      FatalErrorInFunction << "Wrong Ys size" << abort(FatalError);

    if (Tg < 10.0 || Ts < 10.0)
      FatalErrorInFunction << "Too low temperature: \n Tg = "<< Tg << nl << "Ts = " << Ts << abort(FatalError);

    if (rhog < 1e-8)
      FatalErrorInFunction << "Too low density" << abort(FatalError);

    if (p < 10)
      FatalErrorInFunction << "Too low pressure" << abort(FatalError);

    if (ini)
    {
        initializeSurface();
        return 0.0;
    }

    // Calculate reaction rates
    scalarField c0(nSpecie_);
    if (Ts > Treact_)
    {
        rhogs_ = rhog*Tg/Ts;
        for (label i=0; i<nGasSpecie_; i++)
        {
            c_[i] = rhogs_*Yg[i]/gasSpecieThermo_[i].W();
            c0[i] = c_[i];
        }
        for (label i=0; i<nSolidSpecie_; i++)
        {
            c_[i + nGasSpecie_] = Ys[i]*solidSpecieThermo_[i].sden()/solidSpecieThermo_[i].size();
            c0[i + nGasSpecie_] = c_[i + nGasSpecie_];
        }

        // Initialise time progress
        scalar timeLeft = deltaT;

        if (instantaneousRates_)
        {
            omega(p, Ts, c_, celli, c0);

            for (label i=0; i<nGasSpecie_; i++)
            {
                RR[i] = c0[i]*gasSpecieThermo_[i].W();
            }
        }
        else
        {
            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(p, Ts, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            for (label i=0; i<nGasSpecie_; i++)
            {
                RR[i] = (c_[i] - c0[i])*gasSpecieThermo_[i].W()/deltaT;
            }
            for (label i=0; i<nSolidSpecie_; i++)
            {
                //RR[i+nGasSpecie_] = (c_[i+nGasSpecie_] - c0[i+nGasSpecie_])/deltaT
                //        /solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
                Ys[i] = c_[i+nGasSpecie_]/solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
            }
        }
    }
    else
    {
        for (label i=0; i<nGasSpecie_; i++)
        {
            RR[i] = 0.0;
        }
    }

    // Calculate and return reaction heat
    scalar Qdot = 0.0;
    forAll(Yg, i)
    {
        const scalar hi = gasSpecieThermo_[i].Hf();
        Qdot -= hi*RR[i];
    }
    //forAll(Ys, i)
    //{
    //    const scalar hi = solidSpecieThermo_[i].Hf();
    //    Qdot -= hi*RR[i+nGasSpecie_];
    //}
    return Qdot;
}

// ************************************************************************* //
