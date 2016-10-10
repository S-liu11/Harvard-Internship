###########################################################################
#   PURPOSE:
#       INITIATE THE FIRST ROUND OF STEPS
###########################################################################
#   INTERFACE FOR VARIABLES & FUNCTIONS & EXTERNAL FILES(SCRIPT+SUBROUTINE)
#       VARIABLES:
#           SweepRounds, CHECKPOINTS(VARIABLES), NEIGENFREQ,
#       
#       
###########################################################################
#CREATE STEP: GENERAL STATIC-STRECHING
mdb.models['Model-1'].StaticStep(maxNumInc=1000, initialInc=0.001, minInc=1e-6,maxInc=0.1, 
    name= 'Scale-Stretching', nlgeom=ON, previous='Initial')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'U','RF'))

#CREATE FIELDOUTPUT TO CRACK THE HEIGHT OF THE TIP
mdb.models['Model-1'].FieldOutputRequest(createStepName='Scale-Stretching', 
    name='tip_height', rebar=EXCLUDE, region=
    mdb.models['Model-1'].rootAssembly.instances['Part-RVE-1'].sets['tip'], 
    sectionPoints=DEFAULT, variables=('U', ))

NSTRAIN=0.2
NSTRAIN1=NSTRAIN
NSTRAIN2=0.

# Externally Applied Strain Through the Virtual Points (X-DIR) =================
# (engineering strain) =========================================================
# Externally Applied Strain Through the Virtual Points (Y-DIR) =================
# (engineering strain) =========================================================
#   # for the real structure
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', fixed=OFF
    , localCsys=None, name='BC-CornerPinned-real', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-RVE-1'].sets['Set-CornerBL']
    , u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

#   #   # H11, H21
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, 
    createStepName='Scale-Stretching', distributionType=UNIFORM, fieldName=''
    , fixed=OFF, localCsys=None, name='BC-VNReal-H1', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-RealH1'].sets['Set-Point']
    , u1=NSTRAIN1, u2=0, ur3=UNSET)
#   #   # H12, H22
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, 
    createStepName='Scale-Stretching', distributionType=UNIFORM, fieldName=''
    , fixed=OFF, localCsys=None, name='BC-VNReal-H2', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-RealH2'].sets['Set-Point']
    , u1=UNSET, u2=UNSET, ur3=UNSET)

#CREATE STRETCH JOB
Jobname = 'scale_stretching_U3PBC' + str(lencraratio)+'_'+ str(angl)+'_rotate_'+str(dirc)+'_'+str(pot) # Note the power of 1e4 determins the 4 in the brackets
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
#    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    memory=16, memoryUnits=GIGA_BYTES, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=Jobname, nodalOutputPrecision=SINGLE, 
    numCpus=4, numDomains=4, queue=None, scratch='', type=ANALYSIS, 
    userSubroutine='', waitHours=0, waitMinutes=0)
mdb.saveAs(pathName=Jobname)
# LET ME KNOW IF YOU ARE FINISHED ==============================================
print 'Scale_Periodicity_Streching_JOBS.py Accomplished!!!'

#=========================================================================
# 2. SUBMIT THE JOB WITH THE MODIFIED INPUT FILE
#=========================================================================   
mdb.jobs[Jobname].submit(consistencyChecking=OFF)
    # mdb.jobs[CurJobName].submit(consistencyChecking=OFF, datacheckJob=True)
mdb.jobs[Jobname].waitForCompletion()
time.sleep(5)
