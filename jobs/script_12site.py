import os
seeds= ["1"]
#seeds = [1,2,3,4,5,6,7,8,9,10]
#seeds= ["161","162","163","164","165"]
#seeds = [341,342,343,344,345,346,347,348,349,350]
#seeds = [351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370]
#szs=["10", "11","12"]
#
#seeds = [491,492,493,494,495,496,497,498,499,500]
#szs=["-12","-11","-10","-9","-8","-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7","8","9","10","11","12"]
szs = ["0"]
#jbqs=["-5","-4:","-3","-2","-1","-0.5","-0.3","-0.2","-0.1","-0.05","0","0.05","0.1","0.2","0.3","0.5","1","2","3","4","5"]
#jbqs = ["-1","-0.5","-0.2","-0.1","-0.05","0","0.05","0.1","0.2","0.5",1]
jbqs = ["0"]
js=[]
for jbq in jbqs: js.append(["1","1",jbq,jbq])

for j in js:
        jstrong=j[0]
        jweak=j[1]
        jbqstrong=j[2]
        jbqweak=j[3]
        for sz in szs:
                for seed in seeds: #Iterates over objects
                        # Write a input file
                        infname ="in_12site_seed_no_"+str(seed)+"_sz_"+sz+"_jstrong_"+jstrong+"_jweak_"+jweak+"_jbqstrong_"+jbqstrong+"_jbqweak_"+jbqweak+".txt"
                        print infname
                        f = open(infname, 'w')
                        f.write("seed="+str(seed)+"\n")
                        f.write("hamiltonian=spin_model\n")
                        f.write("j_strong="+jstrong+"\n")
                        f.write("j_weak="+jweak+"\n")
                        f.write("jbq_strong="+jbqstrong+"\n")
                        f.write("jbq_weak="+jbqweak+"\n")
                        f.write("spin=1.0\n")
                        f.write("sz="+sz+"\n")
                        f.write("model_strong_triangles=[[0,10,8],[1,5,9],[2,6,4],[3,7,11]]\n")
                        f.write("model_weak_triangles=[[0,1,2],[3,4,5],[6,7,8],[9,10,11]]\n")
                        f.write("[ED]\n")
                        f.write("other_analyses=qddm_spin1\n")
                        f.close()

                        # Write a job file
                        outdirname="/gpfs/home/hchanglani/gittests/FiniteT/jobs/"
                        jobfname ="job_12site_seed_no_"+str(seed)+"_sz_"+sz+"_jstrong_"+jstrong+"_jweak_"+jweak+"_jbqstrong_"+jbqstrong+"_jbqweak_"+jbqweak+".txt"
                        outfname = outdirname+"Hilbert_Space_Avgout_12site_seed_no_"+str(seed)+"_sz_"+sz+"_jstrong_"+jstrong+"_jweak_"+jweak+"_jbqstrong_"+jbqstrong+"_jbqweak_"+jbqweak+".txt"
                        j = open(jobfname, 'w')
                        j.write("#!/bin/bash -l\n")
                        j.write("#Batch Queue Script\n")
                        j.write("#SBATCH --time=4:00:00\n")
                        j.write("#SBATCH --nodes=1\n")
                        j.write("#SBATCH --ntasks-per-node=1\n")
                        j.write("#SBATCH --cpus-per-task=24\n")
                        j.write("#SBATCH -p genacc_q\n")
                        #j.write("#SBATCH -p changlani_q\n")
			j.write("export OMP_NUM_THREADS=24\n")
                        j.write("/gpfs/home/hchanglani/gittests/FiniteT/ed "+infname+" > "+outfname)
                        j.close()
                        os.system("chmod 755 "+jobfname)
                        os.system("sbatch "+jobfname)






	
