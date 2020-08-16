Info << "Reading catalytic walls properties\n" << endl;
IOdictionary catDict
(
    IOobject
    (
        "catalyticWallsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
word catalystName(catDict.lookup("catalystName"));
List<word> catWalls(catDict.lookup("catalyticWalls"));

// Checking boundary conditions
forAll(catWalls, wi)
{
    label catWallID = mesh.boundary().findPatchID(catWalls[wi]);
    forAll(Y, i)
    {
        if (!isA<catalyticWallSpeciesFluxFvPatchScalarField>(Y[i].boundaryField()[catWallID]))
            FatalError << "Boundary field of " << Y[i].name()
                    << " on catalytic wall patch " << catWalls[wi]
                    << " is not a catalyticWallSpeciesFluxFvPatchScalarField."
                    << exit(FatalError);
    }
}

// Write required thermo fields for catalyst 'phase'
{
    forAll(Y,i)
    {
        volScalarField Ysi
        (
            IOobject
            (
                IOobject::groupName(IOobject::member(Y[i].name()), catalystName),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Y[i],
            calculatedFvPatchScalarField::typeName
        );
        Ysi.write();
    }

    volScalarField Tsi
    (
        IOobject
        (
            IOobject::groupName(IOobject::member(thermo.T().name()), catalystName),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.T(),
        calculatedFvPatchScalarField::typeName
    );
    Tsi.write();
}

Info<< "Creating catalyst thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> catThermo(rhoReactionThermo::New(mesh, catalystName));
rhoReactionThermo& cThermo = catThermo();

autoPtr<basicGSChemistryModel> catReaction
(
    basicGSChemistryModel::New(cThermo)
);

forAll(mesh.C(), celli)
{
    cThermo.setCellReacting(celli, false);
    thermo.setCellReacting(celli, true);
}

catReaction->disableSurfaceSpeciesEqn();
