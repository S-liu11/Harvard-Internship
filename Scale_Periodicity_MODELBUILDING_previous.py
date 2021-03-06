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
crack1=(1-lengthcrackratio)/2*BR+(1+lengthcrackratio)/2*TL
crack2=(1+lengthcrackratio)/2*BR+(1-lengthcrackratio)/2*TL
crack3=(1+lengthcrackratio)/2*BR+(1+lengthcrackratio)/2*TL
crack4=(1+lengthcrackratio-0.01)/2*BR+(1+lengthcrackratio-0.01)/2*TL
crack5=(1-lengthcrackratio)/2*BR+(1-lengthcrackratio)/2*TL #this is used to partition face
crack6=(1-lengthcrackratio)/2*BR+TL
crack7=(1+lengthcrackratio)/2*BR+TL
crack8=(1+lengthcrackratio)/2*TL+BR
crack9=(1-lengthcrackratio)/2*TL+BR
crack10=(1+lengthcrackratio)/2*BR
crack11=(1-lengthcrackratio)/2*BR
crack12=(1-lengthcrackratio)/2*TL
crack13=(1+lengthcrackratio)/2*TL
# generate the crack edges
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack1, point2=crack3)   
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack3, point2=crack2) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack2, point2=crack4) 
mdb.models['Model-1'].sketches['__profile__'].Line(point1=crack4, point2=crack1) 

#build the part 
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1UC', type= DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1UC'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

#get a datum point on the position of crack5, this is used to partition face
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack5[0], crack5[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack6[0], crack6[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack7[0], crack7[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack8[0], crack8[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack9[0], crack9[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack10[0], crack10[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack11[0], crack11[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack12[0], crack12[1], 0.0))
mdb.models['Model-1'].parts['Part-1UC'].DatumPointByCoordinate(coords=(crack13[0], crack13[1], 0.0))

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
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack1[0], crack1[1], 
    0.0)), point2=d[2], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack2[0], crack2[1], 
    0.0)), point2=d[2], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack1[0], crack1[1], 
    0.0)), point2=d[3], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack1[0], crack1[1], 
    0.0)), point2=d[10], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack3[0], crack3[1], 
    0.0)), point2=d[4], 
    faces=f)
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack3[0], crack3[1], 
    0.0)), point2=d[5], 
    faces=f.findAt(coordinates=((crack3[0]+crack8[0])/2, (crack3[1]+crack8[1])/2, 
    0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack2[0], crack2[1], 
    0.0)), point2=d[6], 
    faces=f.findAt(coordinates=((crack2[0]+crack9[0])/2, (crack2[1]+crack9[1])/2, 
    0.0)))
p.PartitionFaceByShortestPath(point1=v2.findAt(coordinates=(crack2[0], crack2[1], 
    0.0)), point2=d[7], 
    faces=f.findAt(coordinates=((crack2[0]+crack10[0])/2, (crack2[1]+crack10[1])/2, 
    0.0)))
p.PartitionFaceByShortestPath(point1=d[2], point2=d[8], 
    faces=f.findAt(coordinates=((crack5[0]+crack11[0])/2, (crack5[1]+crack11[1])/2, 
    0.0)))
p.PartitionFaceByShortestPath(point1=d[2], point2=d[9], 
    faces=f.findAt(coordinates=((crack5[0]+crack12[0])/2, (crack5[1]+crack12[1])/2, 
    0.0)))


# mesh
SEEDNUMBER=50
SEEDSIZE1=int(SEEDNUMBER*lengthcrackratio)#SIZE OF THE LARGEST MESH
SEEDSIZE2=int(SEEDNUMBER*(1-lengthcrackratio)/2)
ratio=1#RATIO OF THE SIZE OF THE LARGEST MESH TO THE SMALLEST ONES
#p.seedPart(size=latvec1[0]/SEEDNUMBER, deviationFactor=0.1, minSizeFactor=0.1)
	#seed bias method: DOUBLE, to refine the mesh around the crack corner
	#seed the edges by fixed constraint to guarantee a perfect PBC
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+TL[0])/2,(BL[1]+TL[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+BR[0])/2,(BL[1]+BR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TL[0]+TR[0])/2,(TL[1]+TR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
    constraint=FIXED, endEdges=
    mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TR[0]+BR[0])/2,(TR[1]+BR[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TL[0]+crack6[0])/2,(TL[1]+crack6[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TL[0]+crack13[0])/2,(TL[1]+crack13[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TR[0]+crack7[0])/2,(TR[1]+crack7[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((TR[0]+crack8[0])/2,(TR[1]+crack8[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BR[0]+crack9[0])/2,(BR[1]+crack9[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BR[0]+crack10[0])/2,(BR[1]+crack10[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FINER, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+crack11[0])/2,(BL[1]+crack11[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FIXED, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((BL[0]+crack12[0])/2,(BL[1]+crack12[1])/2,0),)), ratio=ratio, number=SEEDSIZE2)

	#seed the crack and partition edges linked to the crack by finer constraint
	#!WARNING must be FINER, or meshing would be failed!
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FINER, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack2[0]+crack3[0])/2,(crack2[1]+crack3[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FINER, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack2[0]+crack4[0])/2,(crack2[1]+crack4[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FINER, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack1[0]+crack3[0])/2,(crack1[1]+crack3[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
mdb.models['Model-1'].parts['Part-1UC'].seedEdgeByBias(biasMethod=DOUBLE, 
	constraint=FINER, endEdges=
	mdb.models['Model-1'].parts['Part-1UC'].edges.findAt((((crack1[0]+crack4[0])/2,(crack1[1]+crack4[1])/2,0),)), ratio=ratio, number=SEEDSIZE1)
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
    QUAD, regions= mdb.models['Model-1'].parts['Part-1UC'].faces, technique=FREE)
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
for node in mdb.models['Model-1'].parts['Part-RVE'].nodes:
	if (node.coordinates[0]-crack4[0])**2+(node.coordinates[1]-crack4[1])**2 < 1e-10:
		break

# create a set for the tip
mdb.models['Model-1'].parts['Part-RVE'].Set(nodes= 
    (mdb.models['Model-1'].parts['Part-RVE'].nodes[node.label-1:node.label],), name='tip')

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

