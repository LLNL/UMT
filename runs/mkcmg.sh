#!/bin/bash
px=2
py=2
pz=1
let pxm1=$px-1
let pym1=$py-1
let pzm1=$pz-1

echo '#processor block decomposition'
echo 'sms('${px},${py},${pz}')'
echo '#Always specify blocks in block base numbering'
echo 'blk(on,0:'$pxm1',0:'$pym1',0:'$pzm1')'
echo
echo '# tag boundary faces'
echo 'tag("xMinFaces",face,(0:0,0:'$py',0:'$pz'))'
echo 'tag("xMaxFaces",face,('$px':'$px',0:'$py',0:'$pz'))'
echo 'tag("yMinFaces",face,(0:'$px',0:0,0:'$pz'))'
echo 'tag("yMaxFaces",face,(0:'$px','$py':'$py',0:'$pz'))'
echo 'tag("zMinFaces",face,(0:'$px',0:'$py',0:0))'
echo 'tag("zMaxFaces",face,(0:'$px',0:'$py','$pz':'$pz'))'

echo
echo '# define number of zones in each axis'
echo 'numzones('$1','$1','$1')'
echo

echo '#Hex subdivisions'
echo 'sub(10%,0:'$pxm1', 0:'$pym1', 0:'$pzm1',(7,0,0,0)) #7 hex'
echo 'seed(10)'
echo
