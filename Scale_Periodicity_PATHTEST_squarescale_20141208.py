###########################################################################
#   SCALE_PERIODICITY_PATHTEST
#   by SHUANG, version: 20141011  
###########################################################################
#   PURPOSE:
#       DO THE PATH TEST FOR EACH MODEL TO CHECK FOR PATH INDEPENDENCE
###########################################################################


for i in range(0,11):
#CREATE STEP: GENERAL STATIC-STRECHING , WITH PARAMETER
	mdb.models['Model-1'].StaticStep(maxNumInc=1000, initialInc=0.001, minInc=1e-10,maxInc=0.1, 
        name= 'SS_0', nlgeom=ON, previous='Initial')
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial'
        , distributionType=UNIFORM, fieldName='', fixed=OFF
        , localCsys=None, name='BC-CornerPinned-real', region=
        mdb.models['Model-1'].rootAssembly.instances['Part-RVE-1'].sets['Set-CornerBL']
        , u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, 
        createStepName='SS_0', distributionType=UNIFORM, fieldName=''
        , fixed=OFF, localCsys=None, name='BC-VNReal-H1', region=
        mdb.models['Model-1'].rootAssembly.instances['Part-RealH1'].sets['Set-Point']
        , u1=UNSET, u2=UNSET, ur3=UNSET)
	mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, 
        createStepName='SS_0', distributionType=UNIFORM, fieldName=''
        , fixed=OFF, localCsys=None, name='BC-VNReal-H2', region=
        mdb.models['Model-1'].rootAssembly.instances['Part-RealH2'].sets['Set-Point']
        , u1=0.0, u2=0.025*i, ur3=UNSET)
	
#CREATE SHEARING STEP
	for j in range(1,11):
		NSTRAIN2=0.025*j
		pathname='SS_'+str(j)
		previouspathname='SS_'+str(j-1)
		mdb.models['Model-1'].StaticStep(name=pathname, 
            previous=previouspathname, maxNumInc=1000, initialInc=0.001, 
            minInc=1e-10, maxInc=0.1)
		mdb.models['Model-1'].boundaryConditions['BC-VNReal-H1'].setValuesInStep(
	        stepName=pathname, u2=NSTRAIN2)

#CREATE STRETCH JOB
	Jobname = 'scale_stretching_U3PBC_squarescale_tnsh4_' + str(lencraratio)+'_'+ str(perc)+'_PATH_'+str(i)
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
#    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        memory=16, memoryUnits=GIGA_BYTES, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=Jobname, nodalOutputPrecision=SINGLE, 
        numCpus=4, numDomains=4, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)
	mdb.saveAs(pathName=Jobname)
# LET ME KNOW IF YOU ARE FINISHED ==============================================

#=========================================================================
# 2. SUBMIT THE JOB WITH THE MODIFIED INPUT FILE
#=========================================================================   
	mdb.jobs[Jobname].submit(consistencyChecking=OFF)
    # mdb.jobs[CurJobName].submit(consistencyChecking=OFF, datacheckJob=True)
	mdb.jobs[Jobname].waitForCompletion()
	time.sleep(2)

	for k in range(0,11):
		stepname='SS_'+str(k)
		del mdb.models['Model-1'].steps[stepname]

