/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

Class
    Foam::functionObjects::writeUCylindricalCoordinates

Description

    Example of function object specification:
    \verbatim
    writeUCylindricalCoordinates
    {
        type        writeUCylindricalCoordinates;
        libs        ("libcatchyFunctionObjects.so");
        axis        2;
        U           U.gas;
        ...
    }
    \endverbatim

Usage
    \table
        Property  | Description                 | Required    | Default value
        type      | type name: writeUCylindricalCoordinates | yes       |
        axis      | the axis of the cylinder | yes       |
        U         | list of velocity field names | yes |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    writeUCylindricalCoordinates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeUCylindricalCoordinates_H
#define functionObjects_writeUCylindricalCoordinates_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeUCylindricalCoordinates Declaration
\*---------------------------------------------------------------------------*/

class writeUCylindricalCoordinates
:
    public fvMeshFunctionObject
{
    //- Private member variables
    const label axDir;
    const volScalarField rad;
    const volScalarField theta;
    const wordList Unames_;

public:

    //- Runtime type information
    TypeName("writeUCylindricalCoordinates");


    // Constructors

        //- Construct from Time and dictionary
        writeUCylindricalCoordinates
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeUCylindricalCoordinates(const writeUCylindricalCoordinates&) = delete;


    //- Destructor
    virtual ~writeUCylindricalCoordinates();


    // Member Functions

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeUCylindricalCoordinates&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
