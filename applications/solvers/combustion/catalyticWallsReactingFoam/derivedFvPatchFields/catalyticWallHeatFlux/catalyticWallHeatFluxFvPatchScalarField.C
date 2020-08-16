/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------


\*---------------------------------------------------------------------------*/

#include "catalyticWallHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicGSChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyticWallHeatFluxFvPatchScalarField::
catalyticWallHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    catalystName_("catalyst")
{}


Foam::catalyticWallHeatFluxFvPatchScalarField::
catalyticWallHeatFluxFvPatchScalarField
(
    const catalyticWallHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    catalystName_(ptf.catalystName_)
{}


Foam::catalyticWallHeatFluxFvPatchScalarField::
catalyticWallHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    catalystName_(dict.lookupOrDefault<word>("catalystName", "catalyst"))
{}


Foam::catalyticWallHeatFluxFvPatchScalarField::
catalyticWallHeatFluxFvPatchScalarField
(
    const catalyticWallHeatFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    catalystName_(tppsf.catalystName_)
{}


Foam::catalyticWallHeatFluxFvPatchScalarField::
catalyticWallHeatFluxFvPatchScalarField
(
    const catalyticWallHeatFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    catalystName_(tppsf.catalystName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::catalyticWallHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const rhoReactionThermo& thermo = db().lookupObject<rhoReactionThermo>
    (
        IOobject::groupName(basicThermo::dictName, IOobject::group(internalField().name()))
    );
    const basicGSChemistryModel& reaction = db().lookupObject<basicGSChemistryModel>
    (
        IOobject::groupName("chemistryProperties",catalystName_)
    );

    const scalarField Qdotw = reaction.catalyticWallHeatFlux(patchi); // [J/(m2.s)]
    const scalarField kappaw = thermo.kappa(patchi); // [J/(m.s.K)]

    // kappaw*Tgradient = Qdotw
    // [J/(m.s.K)]*[K/m] = [J/(m2.s)]
    gradient() = Qdotw/kappaw; // [K/m]

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::catalyticWallHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);

    os.writeKeyword("catalystName") << catalystName_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        catalyticWallHeatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
