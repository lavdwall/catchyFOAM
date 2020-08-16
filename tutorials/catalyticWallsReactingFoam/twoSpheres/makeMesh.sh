#!/bin/bash -l

rm -rf 0
blockMesh
transformPoints -translate '(-6.0 -2.5 -0.05)'
transformPoints -rollPitchYaw '(0 180 -90)'
mkdir translatedMesh
cp -r constant system translatedMesh
sed -i "s/inertWall/inertWall1/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/reactingWall/reactingWall1/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/axis/axis1/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/inlet/inlet1/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/front/blabla/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/back/front1/g" translatedMesh/constant/polyMesh/boundary
sed -i "s/blabla/back1/g" translatedMesh/constant/polyMesh/boundary


blockMesh
sed -i "s/inertWall/inertWall2/g" constant/polyMesh/boundary
sed -i "s/reactingWall/reactingWall2/g" constant/polyMesh/boundary
sed -i "s/axis/axis2/g" constant/polyMesh/boundary
sed -i "s/inlet/inlet2/g" constant/polyMesh/boundary
sed -i "s/outlet/inlet/g" constant/polyMesh/boundary
sed -i "s/back/back2/g" constant/polyMesh/boundary
sed -i "s/front/front2/g" constant/polyMesh/boundary
transformPoints -translate '(-6.0 -2.5 -0.05)'
transformPoints -rollPitchYaw '(0 0 90)'
mergeMeshes . translatedMesh -overwrite
stitchMesh inlet1 inlet2 -perfect -overwrite
rm -rf 0
rm -r translatedMesh
transformPoints -translate '(-2.5 0 -0.05)'
createPatch -overwrite -dict system/createPatchDict.1
extrudeMesh
collapseEdges -overwrite
transformPoints -scale '(0.001 0.001 0.001)'
createPatch -overwrite -dict system/createPatchDict.2
checkMesh
rm -rf 0
