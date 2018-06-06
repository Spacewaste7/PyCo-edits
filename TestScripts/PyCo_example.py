from PyContact.core.Scripting import PyContactJob, JobConfig

job = PyContactJob("/home/bdv1/Desktop/rpn11_ubq.psf", "/home/bdv1/Desktop/rpn11_ubq.dcd", "attemptPdbDcd1", JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid UBQ", "segid RN11"))

job.runJob(4)

job.writeSessionToFile()
                         
