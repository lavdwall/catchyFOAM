/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "writeUCylindricalCoordinates.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeUCylindricalCoordinates, 0);
    addToRunTimeSelectionTable(functionObject, writeUCylindricalCoordinates, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeUCylindricalCoordinates::writeUCylindricalCoordinates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    axDir(readLabel(dict.lookup("axis"))),
    rad
    (
	    IOobject
	    (
	        "radial_coordinate",
	        runTime.timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    sqrt(sqr(mesh_.C().component((axDir+1)%3))+sqr(mesh_.C().component((axDir+2)%3)))
    ),
    theta
    (
	    IOobject
	    (
	        "tangential_coordinate",
	        runTime.timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    atan2(mesh_.C().component((axDir+1)%3),mesh_.C().component((axDir+2)%3))
    ),
    Unames_(dict.lookup("U"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeUCylindricalCoordinates::~writeUCylindricalCoordinates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeUCylindricalCoordinates::execute()
{
    return true;
}


bool Foam::functionObjects::writeUCylindricalCoordinates::write()
{
    Log << "    Writing cylindrical coordinates"
        << " to " << time_.timeName() << endl;

    rad.write();
    theta.write();

    forAll(Unames_, i)
    {
        volVectorField U = mesh_.lookupObject<volVectorField>(Unames_[i]);

        volScalarField radU
        (
            IOobject
            (
                Unames_[i]+".rad",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (U.component((axDir+1)%3)*mesh_.C().component((axDir+1)%3)+U.component((axDir+2)%3)*mesh_.C().component((axDir+2)%3))/rad
        );
        radU.write();
        volScalarField tanU
        (
            IOobject
            (
                Unames_[i]+".tan",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           (-U.component((axDir+1)%3)*mesh_.C().component((axDir+2)%3)+U.component((axDir+2)%3)*mesh_.C().component((axDir+1)%3))/rad
        );
        tanU.write();
    }

    return true;
}


// ************************************************************************* //
