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

#include "StandardGSChemistryModel.H"
#include "multiComponentGSMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"


// * * * * * * * * * * *  Private member functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
initializeSurface()
{
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
    initialized_ = true;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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
    const scalar rhoi = dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
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

    multiplyAV_ = 0.0;
    multiplyPor_ = 0.0;

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

    multiplyAV_ = 1.0;
    multiplyPor_ = 1.0;

    for (label i=0; i<this->nSolidSpecie_; i++)
    {
        Ys_[i][celli] = c_[i+nGasSpecie_]/solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::StandardGSChemistryModel
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

    nGasSpecie_(Y_.size()),
    nSolidSpecie_(Ys_.size()),
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
    RRg_(nGasSpecie_),
    RRs_(nSolidSpecie_),
    c_(nSpecie_+2),
    dcdt_(nSpecie_+2),
    initialized_(false),
    catarea_(this->lookupOrDefault("catalyticSurfaceAreaPerVolume", 1.0)),
    multiplyAV_(1.0),
    fpor_(this->lookupOrDefault("catalystPorosity", 0.0)),
    por_(fpor_),
    multiplyPor_(1.0),
    RRgwall_(nGasSpecie_),
    RRswall_(nSolidSpecie_)
{
    // Create the fields for the chemistry sources for gas species
    forAll(RRg_, fieldi)
    {
        RRg_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
        RRgwall_.set
        (
            fieldi,
            new volScalarField::Boundary
            (
                this->mesh().boundary(),
                RRg_[fieldi],
                calculatedFvPatchScalarField::typeName
            )
        );
        forAll(RRgwall_[fieldi], facei)
        {
            RRgwall_[fieldi][facei] = 0.0;
        }
    }

    // Create the fields for the chemistry sources for coverages
    forAll(RRs_, fieldi)
    {
        RRs_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Ys_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimless/dimTime, 0)
            )
        );
        RRswall_.set
        (
            fieldi,
            new volScalarField::Boundary
            (
                this->mesh().boundary(),
                RRs_[fieldi],
                calculatedFvPatchScalarField::typeName
            )
        );
        forAll(RRswall_[fieldi], facei)
        {
            RRswall_[fieldi][facei] = 0.0;
        }
    }

    Info<< "StandardGSChemistryModel: " << nl
        << indent << "Number of gas species = " << nGasSpecie_ << nl
        << indent << "Number of solid species = " << nSolidSpecie_ << nl
        << indent << "Number of gas reactions = " << nGasReaction_ << nl
        << indent << "Number of solid reactions = " << nSolidReaction_ << nl
        << indent << "Catalyst porosity = " << fpor_ << nl
        << indent << "Catalyst surface area to volume = " << catarea_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
~StandardGSChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omega
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

    forAll(Y_, i)
    {
        dcdt[i] *= catarea_*multiplyAV_;
    }

    if (por_*multiplyPor_ > 0.0)
    {
        scalarField dcdtg(nGasSpecie_, 0.0);

        forAll(gasReactions_, i)
        {
            const Reaction<GThermoType>& R = gasReactions_[i];

            R.omega(p, T, c, li, dcdtg);
        }

        forAll(dcdtg, i)
        {
            dcdt[i] += dcdtg[i]*por_*multiplyPor_;
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omegaGI
(
    const label index,
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
    const Reaction<GThermoType>& R = gasReactions_[index];
    scalar w = R.omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omegaSI
(
    const label index,
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
    const SurfaceReaction<SThermoType>& R = solidReactions_[index];
    scalar w = R.omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::derivatives
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

    // Constant pressure
    // dT/dt = ...
    /*scalar rho = 0;
    scalar cSum = 0;
    for (label i = 0; i < nGasSpecie_; i++)
    {
        const scalar W = gasSpecieThermo_[i].W();
        cSum += c_[i];
        rho += W*c_[i];
    }
    scalar cp = 0;
    for (label i=0; i<nGasSpecie_; i++)
    {
        cp += c_[i]*gasSpecieThermo_[i].cp(p, T);
    }
    cp /= rho;

    scalar dT = 0;
    for (label i = 0; i < nGasSpecie_; i++)
    {
        const scalar hi = gasSpecieThermo_[i].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;*/

    // dp/dt = ...
    //dcdt[nSpecie_ + 1] = 0;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::jacobian
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
    forAll(Y_, j)
    {
        forAll(c_, i)
        {
            J(j,i) *= catarea_*multiplyAV_;
        }
        dcdt[j] *= catarea_*multiplyAV_;
    }
    if (por_*multiplyPor_>0.0)
    {
        scalarSquareMatrix Jg(nGasSpecie_+1, 0.0);
        scalarField dcdtg(nGasSpecie_+1, 0.0);

        forAll(gasReactions_, gri)
        {
            const Reaction<GThermoType>& R = gasReactions_[gri];
            scalar kfwd, kbwd;
            R.dwdc(p, T, c_, li, Jg, dcdtg, omegaI, kfwd, kbwd, false, dummy);
        }

        forAll(Y_, j)
        {
            forAll(Y_, i)
            {
                J(j,i) += Jg(j,i)*por_*multiplyPor_;
            }
            dcdt[j] += dcdtg[j]*por_*multiplyPor_;
            J(j,nSpecie_) += Jg(j,nGasSpecie_)*por_*multiplyPor_;
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::volScalarField>
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::tc() const
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
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = gasSpecieThermo_[i].Hf();
                Qdot[celli] -= hi*RRg_[i][celli];
            }
        }

        forAll(Ys_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = solidSpecieThermo_[i].Hf();
                Qdot[celli] -= hi*RRs_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::calculateRR
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
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    if (!initialized_)
    {
        initializeSurface();
    }

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(T, celli)
    {
        if (!this->thermo().cellReacting(celli)) continue;

        const scalar Ti = T[celli];
        const scalar pi = p[celli];
        const scalar rhoi = dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).gasCellMixture(celli).rho(pi, Ti);

        for (label i=0; i<nGasSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/gasSpecieThermo_[i].W();
        }
        for (label i=0; i<this->nSolidSpecie_; i++)
        {
            const scalar Ysi = Ys_[i][celli];
            c_[i + nGasSpecie_] = Ysi*solidSpecieThermo_[i].sden()/solidSpecieThermo_[i].size();
        }

        omega(pi, Ti, c_, celli, dcdt_);

        for (label i=0; i<nGasSpecie_; i++)
        {
            RRg_[i][celli] = dcdt_[i]*gasSpecieThermo_[i].W();
        }

        for (label i=0; i<nSolidSpecie_; i++)
        {
            RRs_[i][celli] = dcdt_[i+nGasSpecie_]*solidSpecieThermo_[i].W();
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
template<class DeltaTType>
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    ReactionThermo::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    if (!initialized_)
    {
        initializeSurface();
    }

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c0(nSpecie_);

    forAll(T, celli)
    {
        if (!this->thermo().cellReacting(celli)) continue;

        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            scalar pi = p[celli];
            const scalar rhoi = dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
            (
                this->thermo()
            ).gasCellMixture(celli).rho(pi, Ti);

            for (label i=0; i<nGasSpecie_; i++)
            {
                const scalar Yi = Y_[i][celli];
                c_[i] = rhoi*Yi/gasSpecieThermo_[i].W();
                c0[i] = c_[i];
            }
            for (label i=0; i<this->nSolidSpecie_; i++)
            {
                const scalar Ysi = Ys_[i][celli];
                c_[i + nGasSpecie_] = Ysi*solidSpecieThermo_[i].sden()/solidSpecieThermo_[i].size();
                c0[i + nGasSpecie_] = c_[i + nGasSpecie_];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<nGasSpecie_; i++)
            {
                RRg_[i][celli] =
                    (c_[i] - c0[i])*gasSpecieThermo_[i].W()/deltaT[celli];
            }

            for (label i=0; i<nSolidSpecie_; i++)
            {
                if (!this->thermo().composition().active(i+nGasSpecie_))
                {
                    Ys_[i][celli] =
                        c_[i+nGasSpecie_]/solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
                }
                RRs_[i][celli] =
                    (c_[i + nGasSpecie_] - c0[i+nGasSpecie_])/deltaT[celli]
                    /solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
            }
        }
        else
        {
            for (label i=0; i<nGasSpecie_; i++)
            {
                RRg_[i][celli] = 0;
            }
            for (label i=0; i<nSolidSpecie_; i++)
            {
                RRs_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
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
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setGasMassFractions
(
    const PtrList<volScalarField>& Yg
)
{
    forAll(Y_, i)
    {
        forAll(Y_[i], celli)
        {
            Y_[i][celli] = Yg[i][celli];
        }
        Y_[i].correctBoundaryConditions();
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
disableSurfaceSpeciesEqn()
{
    Info << "Disabling solution of surface species equation" << nl
        << indent << "NOTE: this means that values of surface species fractions "
        << "will be assigned after solution of the chemistry" << endl;
    forAll(Ys_, i)
    {
        dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
        (
            this->thermo()
        ).setInactive(i+nGasSpecie_);
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
correctRates
(
    PtrList<volScalarField>& Yg,
    const volScalarField& rhog,
    const volScalarField& Tg,
    const volScalarField& pg
)
{
    if (!this->chemistry_)
    {
        Info << "Chemistry is off" << endl;
        return;
    }
    Info << "Correcting catalytic zone rates" << endl;

    setGasMassFractions(Yg);
    scalarField& T = this->thermo().T().primitiveFieldRef();
    scalarField& p = this->thermo().p().primitiveFieldRef();

    multiplyAV_ = 1.0;
    multiplyPor_ = 1.0;

    forAll(T, celli)
    {
        T[celli] = Tg[celli];
        p[celli] = pg[celli];
    }

    if (this->lookupOrDefault("instantaneousRates", false))
    {
        Info << "Calculating instantaneous rates" << endl;
        this->calculate();
    }
    else
    {
        this->solve(this->mesh().time().deltaTValue());
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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
    if (!this->chemistry_)
    {
        Info << "Wall chemistry is off" << endl;
        return;
    }

    if (!initialized_)
    {
        initializeSurface();
    }

    Info << "Correcting catalytic wall fluxes" << endl;

    scalarField c(nSpecie_, 0.0);
    scalarField c0(nSpecie_, 0.0);

    scalar deltaT = this->mesh().time().deltaTValue();

    forAll(patches, patchii)
    {
        label patchi = this->mesh().boundary().findPatchID(patches[patchii]);
        const fvPatchScalarField& pT = Tg.boundaryField()[patchi];
        forAll(pT, facei)
        {
            label celli = pT.patch().faceCells()[facei];
            multiplyAV_ = pT.patch().magSf()[facei]/this->mesh().V()[celli];
            multiplyPor_ = 0.0;

            scalar Ti = Tg.boundaryField()[patchi][facei];
            scalar pi = pg.boundaryField()[patchi][facei];
            scalar rhoi = rhog.boundaryField()[patchi][facei];
            this->thermo().T().boundaryFieldRef()[patchi][facei] = Ti;
            this->thermo().p().boundaryFieldRef()[patchi][facei] = pi;

            forAll(Yg, i)
            {
                RRgwall_[i][patchi][facei] = 0.0;
                scalar Ygi = Yg[i].boundaryField()[patchi][facei];
                Y_[i].boundaryFieldRef()[patchi][facei] = Ygi;
                c[i] = rhoi*Ygi/gasSpecieThermo_[i].W();
                c0[i] = c[i];
            }

            for (label i=0; i<this->nSolidSpecie_; i++)
            {
                scalar Ysi = Ys_[i].boundaryField()[patchi][facei];
                c[i + nGasSpecie_] = Ysi*solidSpecieThermo_[i].sden()/solidSpecieThermo_[i].size();
                c0[i + nGasSpecie_] = c[i + nGasSpecie_];
            }

            // Initialise time progress
            scalar timeLeft = deltaT;

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<this->nGasSpecie_; i++)
            {
                RRgwall_[i][patchi][facei] =
                        (c[i] - c0[i])*this->gasSpecieThermo_[i].W()/deltaT/multiplyAV_;
            }

            scalar sumCS = 0.0;
            for (label i=1; i<this->nSolidSpecie_; i++)
            {
                Ys_[i].boundaryFieldRef()[patchi][facei] = max(c[i + nGasSpecie_]
                    /solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size(),0.0);
                sumCS += Ys_[i].boundaryFieldRef()[patchi][facei];
            }
            Ys_[0].boundaryFieldRef()[patchi][facei] = 1.0 - sumCS;
            for (label i=0; i<this->nSolidSpecie_; i++)
            {
                RRswall_[i][patchi][facei] =
                    (c[i+nGasSpecie_] - c0[i+nGasSpecie_])/deltaT;
            }
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallSpeciesFlux
(
    const word speciei,
    const label patchi
) const
{
    label i = this->thermo().composition().species()[speciei];
    return RRgwall_[i][patchi];
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallHeatFlux
(
    const label patchi
) const
{
    scalarField tQdot(this->mesh().boundary()[patchi].size(), 0.0);

    if (this->chemistry_)
    {
        forAll(Y_, i)
        {
            forAll(tQdot, facei)
            {
                scalar hf = gasSpecieThermo_[i].Hf();
                tQdot[facei] -= hf*RRgwall_[i][patchi][facei];
            }
        }
        forAll(Ys_, i)
        {
            forAll(tQdot, facei)
            {
                scalar hf = solidSpecieThermo_[i].Hf();
                tQdot[facei] -= hf*RRswall_[i][patchi][facei];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setAlpha
(
    const volScalarField& alpha
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTemperature
(
    const volScalarField& T
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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
    const bool
)
{
    Info << "Select gsParticle chemistry type if you want to use this model" << endl;
    NotImplemented;
    return 0.0;
}

// ************************************************************************* //
