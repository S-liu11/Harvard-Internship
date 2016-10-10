###########################################################################
#   Scale_Periodicity_MODELBUILDING
#   BUILD PART, INSTANCE AND MESH    
#       
#       
#       
###########################################################################
# Lattice Vector
latvec1 = numpy.array([10.,0.])
latvec2 = numpy.array([10.*cos(ang),10.*sin(ang)])
BL = numpy.array([0.,0.])
BR = BL + 1*latvec1
TL = BL + 1*latvec2
TR = BL + 1*latvec1 + 1*latvec2

# sketch for a single RVE
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=300.0)

# generate the edges  
mdb.models['Model-1'].sketches['__profile__'].Line(point1=BL, point2=TL)   
mdb.models['Model-1'].sketches['__profile__'].Line(point1=TL, point2=TR) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=TR, point2=BR) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=BR, point2=BL) 

#generate some feature points to build the crack and mesh
crack1=((1-lengthcrackratio)/2*BR+(1+lengthcrackratio)/2*TL)/2
crack2=((1+lengthcrackratio)/2*BR+(1-lengthcrackratio)/2*TL)/2
crack3=((1+lengthcrackratio)/2*BR+(1+lengthcrackratio)/2*TL)/2
crack4=((1+lengthcrackratio-0.01)/2*BR+(1+lengthcrackratio-0.01)/2*TL)/2
crack5=((1-lengthcrackratio+0.01)/2*BR+(1+lengthcrackratio-0.01)/2*TL)/2 #this is used to partition face
crack6=((1+lengthcrackratio-0.01)/2*BR+(1-lengthcrackratio+0.01)/2*TL)/2
crack7=((1-lengthcrackratio)/2*BR+((1-lengthcrackratio)/2+conratio)*TL)/2
crack8=((1-lengthcrackratio)/2*TL+((1-lengthcrackratio)/2+conratio)*BR)/2
crack9=((1-lengthcrackratio)/2*BR+(1-lengthcrackratio)/2*TL)/2
krack1=crack1+(TL+BR)/2
krack2=crack2+(TL+BR)/2
krack3=crack3+(TL+BR)/2
krack4=crack4+(TL+BR)/2
krack5=crack5+(TL+BR)/2
krack6=crack6+(TL+BR)/2
krack7=crack7+(TL+BR)/2
krack8=crack8+(TL+BR)/2
krack9=crack9+(TL+BR)/2
cen=(TL+BR)/2
# generate the crack edges
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack1, point2=crack3)   
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack3, point2=crack2) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack2, point2=crack8) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack8, point2=crack6) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack6, point2=crack4)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack4, point2=crack5)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack5, point2=crack7)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack7, point2=crack1)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack1, point2=krack3)   
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack3, point2=krack2) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack2, point2=krack8) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack8, point2=krack6) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack6, point2=krack4)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack4, point2=krack5)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack5, point2=krack7)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=krack7, point2=krack1)

#build the part 
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1UC', type= DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1UC'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

#get a datum point on the position of crack5, this is used to partition face
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack9[0], crack9[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(krack9[0], krack9[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(BR[0]/2, BR[1]/2, 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(TL[0]/2, TL[1]/2, 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(TL[0]/2+BR[0], TL[1]/2+BR[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(TL[0]+BR[0]/2, TL[1]+BR[1]/2, 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(cen[0], cen[1], 0.0))
#Part module: Virtual Nodes!!
    # part 2,3: virtual nodes.  Each should have 2 Nodes of 2 DOFs, F(2,2)
    #   Part-F1: Ux->F11,  Uy->F21,  shared by real and imaginary mesh
    #   Part-F2: Ux->F12,  Uy->F22,  shared by real and imaginary mesh
mdb.models['Model-1'].Part(dimensionality=THREE_D, 
    name='Part-Point', type=DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-Point'].ReferencePoint(point=(0.0, 0.0, 0.0))
    # Assign a set to the virtual node
mdb.models['Model-1'].parts['Part-Point'].Set(name='Set-Point', referencePoints=(
    mdb.models['Model-1'].parts['Part-Point'].referencePoints[1], )) 

#partition face to get a better meshing quality
p = mdb.models['Model-1'].parts['Part-1UC']
f = p.faces
v2, e, d = p.vertices, p.edges, p.datums
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack7[0], crack7[1], 
    0.0)), point2=d[2], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack8[0], crack8[1], 
    0.0)), point2=d[2], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack7[0], krack7[1], 
    0.0)), point2=d[3], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack8[0], krack8[1], 
    0.0)), point2=d[3], 
    faces=f)
p.PartitionFaceByShortestPath(point1=d[4], point2=d[8], faces=f)
p.PartitionFaceByShortestPath(point1=d[5], point2=d[8], faces=f)
p.PartitionFaceByShortestPath(point1=d[6], point2=d[8], faces=f)
p.PartitionFaceByShortestPath(point1=d[7], point2=d[8], faces=f.findAt(coordinates=(BR[0]*0.5+TL[0]*0.75, BR[1]*0.5+TL[1]*0.75, 
    0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack9[0], crack9[1], 
    0.0)), point2=v2.findAt(coordinates=(BL[0], BL[1], 0.0)),
    faces=f.findAt(coordinates=((BL[0]+crack9[0])/2, (BL[1]+crack9[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack2[0], crack2[1], 
    0.0)), point2=v2.findAt(coordinates=(BR[0]/2, BR[1]/2, 0.0)),
    faces=f.findAt(coordinates=((BR[0]/2+crack2[0])/2, (BR[1]/2+crack2[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack1[0], crack1[1], 
    0.0)), point2=v2.findAt(coordinates=(TL[0]/2, TL[1]/2, 0.0)),
    faces=f.findAt(coordinates=((TL[0]/2+crack1[0])/2, (TL[1]/2+crack1[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack3[0], crack3[1], 
    0.0)), point2=v2.findAt(coordinates=(cen[0], cen[1], 0.0)),
    faces=f.findAt(coordinates=((cen[0]+crack3[0])/2, (cen[1]+crack3[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack9[0], krack9[1], 
    0.0)), point2=v2.findAt(coordinates=(BL[0]+cen[0], BL[1]+cen[1], 0.0)),
    faces=f.findAt(coordinates=((BL[0]+krack9[0]+cen[0])/2, (BL[1]+krack9[1]+cen[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack2[0], krack2[1], 
    0.0)), point2=v2.findAt(coordinates=(BR[0]/2+cen[0], BR[1]/2+cen[1], 0.0)),
    faces=f.findAt(coordinates=((BR[0]/2+krack2[0]+cen[0])/2, (BR[1]/2+krack2[1]+cen[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack1[0], krack1[1], 
    0.0)), point2=v2.findAt(coordinates=(TL[0]/2+cen[0], TL[1]/2+cen[1], 0.0)),
    faces=f.findAt(coordinates=((TL[0]/2+krack1[0]+cen[0])/2, (TL[1]/2+krack1[1]+cen[1])/2, 0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(krack3[0], krack3[1], 
    0.0)), point2=v2.findAt(coordinates=(cen[0]+cen[0], cen[1]+cen[1], 0.0)),
    faces=f.findAt(coordinates=((cen[0]+krack3[0]+cen[0])/2, (cen[1]+krack3[1]+cen[1])/2, 0.0)))
# mesh
SEEDNUMBER=50
SEEDSIZE1=int(SEEDNUMBER*lengthcrackratio)#SIZE OF THE LARGEST MESH
SEEDSIZE2=int(SEEDNUMBER*(1-lengthcrackratio)/2)
SEEDSIZE3=int(SEEDNUMBER*(lengthcrackratio-conratio))
ratio=1#RATIO OF THE SIZE OF THE LARGEST MESH TO THE SMALLEST ONES
#p.seedPart(size=latvec1[0]/SEEDNUMBER, deviationFactor=0.1, minSizeFactor=0.1)
	#seed bias method: DOUBLE, to refine the mesh around the crack corner
	#seed the edges by fixed constraint to guarantee a perfect PBC
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+TL[0])/2,(BL[1]+TL[1])*0.25,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+BR[0])*0.25,(BL[1]+BR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TL[0]+TR[0])*0.25,(TL[1]+TR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TR[0]+BR[0])/2,(TR[1]+BR[1])*0.25,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+TL[0])/2,(BL[1]+TL[1])*0.75,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+BR[0])*0.75,(BL[1]+BR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TL[0]+TR[0])*0.75,(TL[1]+TR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TR[0]+BR[0])/2,(TR[1]+BR[1])*0.75,0),)), ratio=ratio, number=SEEDSIZE1)

#mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
#	constraint=FINER, endEdges=
#	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack1[0]+crack5[0])/2,(crack1[1]+crack5[1])/2,0),)), maxSize=SEEDSIZE*lengthcrackratio, minSize=SEEDSIZE*lengthcrackratio/ratio)
#mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
#	constraint=FINER, endEdges=
#	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack2[0]+crack5[0])/2,(crack2[1]+crack5[1])/2,0),)), maxSize=SEEDSIZE*lengthcrackratio, minSize=SEEDSIZE*lengthcrackratio/ratio)
#mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
#	constraint=FINER, endEdges=
#	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack1[0]+TL[0])/2,(crack1[1]+TL[1])/2,0),)), maxSize=SEEDSIZE*lengthcrackratio, minSize=SEEDSIZE*lengthcrackratio/ratio)
#mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
#	constraint=FINER, endEdges=
#	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack2[0]+BR[0])/2,(crack2[1]+BR[1])/2,0),)), maxSize=SEEDSIZE*lengthcrackratio, minSize=SEEDSIZE*lengthcrackratio/ratio)


	#mesh control: QUAD, STRUCTURED, S4R
#mdb.models['Model-1'].parts['Part-1UC'].setMeshControls(elemShape=
#    QUAD, regions= mdb.models['Model-1'].parts['Part-1UC'].faces, technique=STRUCTURED) #QUAD
mdb.models['Model-1'].parts['Part-1UC'].setMeshControls(elemShape=
    QUAD, regions= mdb.models['Model-1'].parts['Part-1UC'].faces, technique=STRUCTURED)
#mdb.models['Model-1'].parts['Part-1UC'].setMeshControls(elemShape=
#    TRI, regions= mdb.models['Model-1'].parts['Part-1UC'].faces, technique=FREE)
mdb.models['Model-1'].parts['Part-1UC'].setMeshControls(algorithm=
    DEFAULT, regions= mdb.models['Model-1'].parts['Part-1UC'].faces)
mdb.models['Model-1'].parts['Part-1UC'].setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=STANDARD),), regions=(
   mdb.models['Model-1'].parts['Part-1UC'].faces, ))
#mdb.models['Model-1'].parts['Part-1UC'].setElementType(elemTypes=(ElemType(
#    elemCode=S3, elemLibrary=STANDARD),), regions=(
#    mdb.models['Model-1'].parts['Part-1UC'].faces, ))

# GENERATE MESH FIRST ==========================================================
mdb.models['Model-1'].parts['Part-1UC'].generateMesh()

# creat meshpart ===============================================================
mdb.models['Model-1'].parts['Part-1UC'].PartFromMesh(name='Part-RVE')
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

# Instance it into assembly =================
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, 
    name='Part-RVE-1', part=mdb.models['Model-1'].parts['Part-RVE'])
	
# create set all ===============================================================
mdb.models['Model-1'].parts['Part-RVE'].Set(nodes= 
    mdb.models['Model-1'].parts['Part-RVE'].nodes, name='Set-all')

# find out the tip of the crack. this is used as a criterion of the deformation capability 
for node1 in mdb.models['Model-1'].parts['Part-RVE'].nodes:
	if (node1.coordinates[0]-crack4[0])**2+(node1.coordinates[1]-crack4[1])**2 < 1e-10:
		break

# create a set for the tip
mdb.models['Model-1'].parts['Part-RVE'].Set(nodes= 
    (mdb.models['Model-1'].parts['Part-RVE'].nodes[node1.label-1:node1.label],), name='tip')

for node2 in mdb.models['Model-1'].parts['Part-RVE'].nodes:
	if (node2.coordinates[0]-crack5[0])**2+(node2.coordinates[1]-crack5[1])**2 < 1e-10:
		break

# create a set for the sidetip1
mdb.models['Model-1'].parts['Part-RVE'].Set(nodes= 
    (mdb.models['Model-1'].parts['Part-RVE'].nodes[node2.label-1:node2.label],), name='sidetip1')

for node3 in mdb.models['Model-1'].parts['Part-RVE'].nodes:
	if (node3.coordinates[0]-crack6[0])**2+(node3.coordinates[1]-crack6[1])**2 < 1e-10:
		break

# create a set for the sidetip2
mdb.models['Model-1'].parts['Part-RVE'].Set(nodes= 
    (mdb.models['Model-1'].parts['Part-RVE'].nodes[node3.label-1:node3.label],), name='sidetip2')

#=================
#Material module: EliteDouble32 (light green)
#=================
# density ======================================================================
RHO=1.0E-9;
mdb.models['Model-1'].Material(name='Material-rubber')
mdb.models['Model-1'].materials['Material-rubber'].Density(table=((RHO, ), ))

# Elastic property: Neo-Hookean ================================================
G0=0.2616
NU0=0.4950      # from K0/G0=100

C10=G0/2.0
D1=3.0*(1.0-2.0*NU0) / (G0 *(1.0+2.0*NU0))
 
mdb.models['Model-1'].materials['Material-rubber'].Hyperelastic(
    table=((C10, D1), ), 
    testData=OFF, type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

#=========================================================================
# SECTION
#=========================================================================
# HomogeneousSolidSection ======================================================
mdb.models['Model-1'].HomogeneousShellSection(material='Material-rubber', name=
    'Section-rubber', thickness=0.1)

# SECTION ASSIGNMENT ===========================================================
#mdb.models['Model-1'].parts['Part-lurc'].SectionAssignment(
#    offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
#    mdb.models['Model-1'].parts['Part-lurc'].sets['Set-all']
#    , sectionName='Section-rubber')    # this command is limited to the not-previously-meshed parts
mdb.models['Model-1'].parts['Part-RVE'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    elements=mdb.models['Model-1'].parts['Part-RVE'].elements), 
    sectionName='Section-rubber')

# LET ME KNOW IF YOU ARE FINISHED ==============================================
print 'Material definition Accomplished!!!'




#=================
# 2. Assembly: INSTANCE VIRTUAL NODES INTO ASSEMBLY
#=================
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-RealH1',
    part=mdb.models['Model-1'].parts['Part-Point'])
mdb.models['Model-1'].rootAssembly.instances['Part-RealH1'].translate(
    vector=(0.0, 0.0, 0.0))
#   # FOR THE H12, H22
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-RealH2',
    part=mdb.models['Model-1'].parts['Part-Point'])
mdb.models['Model-1'].rootAssembly.instances['Part-RealH2'].translate(
    vector=(0.0, 0.0, 0.0))
    
    
    
# LET ME KNOW IF YOU ARE FINISHED ==============================================
print 'Scale_Periodicity_MODELBUILDING_20141011.py Accomplished!!!'

