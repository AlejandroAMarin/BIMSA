import os
import sys
import pandas as pd
import numpy as np


rootdir =  os.path.dirname(os.path.realpath(__file__)) + "/../upmem"

resultdir = rootdir + "/profile"

applications = ["wfa"]
def get_data(app_name):
    os.chdir(rootdir)

    chart = open(rootdir + "/../results/" + app_name + "_upmem_output.csv","w")
    chart.write("NumPairs,SeqLength,Error,DPU,TL,BL,BLI,LoadTime,KernelTime,RetriveTime,CPUTime\n")

    #for subdir, dirs, files in os.walk(rootdir+"/"+resultdir):
    for file in os.listdir(resultdir):
        if (not file.endswith("out")):
            continue

        np = file.split("_")[1].replace("n","")
        sl = file.split("_")[2].replace("l","")
        er = file.split("_")[3].replace("e","")
        tl = file.split("_")[4].replace("tl","")
        bl = file.split("_")[5].replace("bl","")
        bli = file.split("_")[6].replace("bli","")
        dpus = file.split("_")[7].replace("dpus","")
        dpus = dpus[:-4]

        avg_load_time = [0]
        avg_kernel_time = [0]
        avg_retrieve_time = [0]
        avg_cpu_time = [0]
        #tmp =  os.path.join(subdir, file)
        tmp =  resultdir + "/" + file
        with open(tmp, "r") as ins:
            for line in ins:
                if(line.find("CPU-DPU (ms):")!=-1):
                    value =  (float(line.split("CPU-DPU (ms): ")[1].split()[0]))
                    avg_load_time.append(value)
                if(line.find("DPU Kernel (ms):")!=-1):
                    value =  (float(line.split("DPU Kernel (ms): ")[1].split()[0]))
                    avg_kernel_time.append(value)
                if(line.find("DPU-CPU Time (ms):")!=-1):
                    value =  (float(line.split("DPU-CPU Time (ms): ")[1].split()[0]))
                    avg_retrieve_time.append(value)
                if(line.find("CPU Version (ms):")!=-1):
                    value =  (float(line.split("CPU Version (ms): ")[1].split()[0]))
                    avg_cpu_time.append(value)

        try:
            chart.write(str(np) + "," + str(sl) + "," + str(er) + "," + str(dpus) + "," + str(tl) + "," + str(bl) + "," + str(bli) + "," + str(max(avg_load_time)) + "," + str(max(avg_kernel_time))+","+str(max(avg_retrieve_time))+"," + str(max(avg_cpu_time))+ "\n")
        except:
            print(tmp)
        
    chart.close()



def main():
    if(len(sys.argv) < 2):
        print ("Usage: python run.py application")
        print ("Applications available: ")
        for value in applications:
            print (value)
        print ("All")
    else:
        cmd = sys.argv[1]
        print ("Application to run is: " + cmd)
        get_data(cmd)

if __name__ == "__main__":
    main()
