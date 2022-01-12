# This is a pmesh mesh generation script.  It can be run with the LLNL parallel mesh generator code 'pmesh'.
__doc__ = """This example produces a crooked pipe mesh, which looks something like:

      rMax-> +------------------------------------------------------+
             |                                                      |
             |                                                      |
             |                                                      |
             |       rBend1-> +---------------------+               |
             |                |                     |               |
             |                |                     |               |
             |                |                     |               |
             |    rInnerEnd-> |     +---------+     |               |
             |                |     |         |     |               |
             |                |     |         |     |               |
    rBend0-> +----------------+     |         |     +---------------+
             |                ^     |         |     ^               |
             |                |     |         |     |               |
             |              zBend0  |         |   zBend1            |
      rMin-> +----------------------+---------+---------------------+
             ^                      ^         ^                     ^
             |                      |         |                     |
            zMin            zInnerStart     zInnerEnd              zMax
"""

## Initialize the Problem Environment ##

import math
Geometry.mode = "rz"

zMin        = 0.
zBend0      = 2.5
zInnerStart = 3.
zInnerEnd   = 4.
zBend1      = 4.5
zMax        = 7.

rMin      = 0.
rBend0    = 0.5
rInnerEnd = 1.
rBend1    = 1.5
rMax      = 2.

## Create the Problem Geometry ##

innerSheet = Sheet.FromPoints( [
    PositionRZ( rMin, zInnerEnd), PositionRZ( rInnerEnd, zInnerEnd ),
    PositionRZ( rInnerEnd, zInnerStart), PositionRZ( rMin, zInnerStart )
] )

middleSheet = Sheet.FromPoints( [
    PositionRZ( rMin, zMin ), PositionRZ( rMin, zInnerStart ),
    PositionRZ( rInnerEnd, zInnerStart ), PositionRZ( rInnerEnd, zInnerEnd ),
    PositionRZ( rMin, zInnerEnd ), PositionRZ( rMin, zMax ),
    PositionRZ( rBend0, zMax ), PositionRZ( rBend0, zBend1 ),
    PositionRZ( rBend1, zBend1 ), PositionRZ( rBend1, zBend0 ),
    PositionRZ( rBend0, zBend0 ), PositionRZ( rBend0, zMin )
] )

outerSheet = Sheet.FromPoints( [
    PositionRZ( rBend0, zMin ), PositionRZ( rBend0, zBend0 ),
    PositionRZ( rBend1, zBend0 ), PositionRZ( rBend1, zBend1 ),
    PositionRZ( rBend0, zBend1 ), PositionRZ( rBend0, zMax ),
    PositionRZ( rMax, zMax ), PositionRZ( rMax, zMin )
] )

globInner = Glob( innerSheet, name="inner", region="thick", color="red", density = 1.59 )
globMiddle = Glob( middleSheet, name="middle", region="thin", color="orange", density = 1.60 )
globOuter = Glob( outerSheet, name="outer", region="thick", color="red", density = 1.61 )

## Annotate the Problem Geometry ##

mfem  =  { "thin":  1,  "thick":  2 }

sourceEdge = globMiddle.borderEdges[-1]

for edge in Geometry.getLinearEdges( r = rMin ):
    edge.addTag( MeshEdgeTag("bottomFaces", "int"), 14 )

for edge in Geometry.getLinearEdges( r = rMax ):
    edge.addTag( MeshEdgeTag("topFaces", "int"), 12 )

for edge in Geometry.getLinearEdges( z = zMax ):
    edge.addTag( MeshEdgeTag("rightFaces", "int"), 13 )

for edge in Geometry.getLinearEdges( z = zMin ):
    if edge != sourceEdge:
        edge.addTag( MeshEdgeTag("leftFaces", "int"), 11 )

sourceEdge.addTag( MeshEdgeTag("sourceFaces", "int"), 10 )

## Webcut the Problem Geometry ##

# NOTE: Webcutting this geometry serves three primary functions:
# 1. Controlling the mesh resolution in the inner and outer regions near the middle region.
# 2. Decomposing the problem geometry into meshable quad faces.
# 3. Introducing virutal edges so that the geometry can be swept.

webcutOffset = 0.1 # Distance of increased resolution in inner and outer

# (1) Webcut the Inner/Outer Near Middle #

innerNearMiddleRegionWire = Wire.FromPoints( [
    PositionRZ( rMin, zInnerStart + webcutOffset ),
    PositionRZ( rInnerEnd - webcutOffset, zInnerStart + webcutOffset ),
    PositionRZ( rInnerEnd - webcutOffset, zInnerEnd - webcutOffset ),
    PositionRZ( rMin, zInnerEnd - webcutOffset )
] )
innerNearMiddleRegionTool = WebcutTool.FromWire( innerNearMiddleRegionWire )
globInner.webcut( innerNearMiddleRegionTool )
innerNearMiddleRegionTool.delete()
innerNearMiddleRegionWire.delete()

outerNearMiddleRegionWire = Wire.FromPoints( [
    PositionRZ( rBend0 + webcutOffset, zMin ),
    PositionRZ( rBend0 + webcutOffset, zBend0 - webcutOffset ),
    PositionRZ( rBend1 + webcutOffset, zBend0 - webcutOffset ),
    PositionRZ( rBend1 + webcutOffset, zBend1 + webcutOffset ),
    PositionRZ( rBend0 + webcutOffset, zBend1 + webcutOffset ),
    PositionRZ( rBend0 + webcutOffset, zMax )
] )
outerNearMiddleRegionTool = WebcutTool.FromWire( outerNearMiddleRegionWire )
globOuter.webcut( outerNearMiddleRegionTool )
outerNearMiddleRegionTool.delete()
outerNearMiddleRegionWire.delete()

# (2) Webcut the Geometry into Quad Faces #

leftOuterBendTool = WebcutTool.Line(
    PositionRZ( rInnerEnd - webcutOffset, zInnerStart + webcutOffset ),
    PositionRZ( rBend1 + webcutOffset, zBend0 - webcutOffset )
)
rightOuterBendTool = WebcutTool.Line(
    PositionRZ( rInnerEnd - webcutOffset, zInnerEnd - webcutOffset ),
    PositionRZ( rBend1 + webcutOffset, zBend1 + webcutOffset )
)

for outerBendTool in [ leftOuterBendTool, rightOuterBendTool ]:
    globInner.webcut( outerBendTool )
    globOuter.webcut( outerBendTool )
    outerBendTool.delete()

leftInnerBendTool = WebcutTool.Line(
    PositionRZ( rMin, zInnerStart ),
    PositionRZ( rBend0 + webcutOffset, zBend0 - webcutOffset )
)
rightInnerBendTool = WebcutTool.Line(
    PositionRZ( rMin, zInnerEnd ),
    PositionRZ( rBend0 + webcutOffset, zBend1 + webcutOffset )
)

for innerBendTool in [ leftInnerBendTool, rightInnerBendTool ]:
    globOuter.webcut( innerBendTool )
    innerBendTool.delete()

propagationVertices = [
    globMiddle.borderVertices[2], globMiddle.borderVertices[3],
    globMiddle.borderVertices[7], globMiddle.borderVertices[10],
    globOuter.vertices[10], globOuter.vertices[12],
]
propagationTools = []

for propagationVertex in propagationVertices:
    propagationEdges = [ e for e in propagationVertex.edges if \
        e.start.position.r == e.end.position.r or \
        e.start.position.z == e.end.position.z ]

    for propagationEdge in propagationEdges:
        edgeStart = propagationEdge.end if propagationEdge.start == propagationVertex \
            else propagationEdge.start
        edgeVector = edgeStart.position.getDirection( propagationVertex.position )

        propagationTool = WebcutTool.Line(
            propagationVertex.position,
            propagationVertex.position + 10.0 * edgeVector,
        )
        propagationTools.append( propagationTool )

for propagationTool in propagationTools:
    for glob in [ globInner, globMiddle, globOuter ]:
        glob.webcut( propagationTool )
    propagationTool.delete()

# (3) Introducing Virtual Edges for Sweep #

virtualCoreTool = WebcutTool.Line(
    PositionRZ( (rMin + rBend0) / 2.0, zMin ),
    PositionRZ( (rMin + rBend0) / 2.0, zMax )
)
globMiddle.webcut( virtualCoreTool )
globInner.webcut( virtualCoreTool )
virtualCoreTool.delete()

## Add Zone Information to the Problem Geometry ##

for edge in Geometry.getLinearEdges( r=0 ) + Geometry.getLinearEdges( z=0 ):
    edge.numZones = math.ceil( 10 * edge.length )

ratioSpacing = ZoneSpacing.Geometric( ratio = 1.47 )
globInner.borderEdges[0].numZones = 10             # inner ratio zoning section
globInner.borderEdges[0].setZoneSpacing( ratioSpacing, globInner.borderVertices[0] )
globOuter.borderEdges[-1].numZones = 10            # outer ratio zoning section
globOuter.borderEdges[-1].setZoneSpacing( ratioSpacing, globOuter.borderVertices[0] )

## Mesh the Geometry and Output Results ##

blocking = Blocking( )
Superblock.createSuperfaces( )
blocking.save( "TopHatMesh.sat" )
mesh = Mesh( )
mesh.saveNURBs( mfem = "crooked_pipe.mesh", materials = mfem, elementSource = 'zone_linear' )
