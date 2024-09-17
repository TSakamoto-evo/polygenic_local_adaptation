import os
import subprocess
import sys
import math
import time

args = sys.argv

para_list = open("parameter_list_full.csv", mode="r")

line_no = 0
for lines in para_list:
    if line_no == int(args[1]):
        line = lines.rstrip("\n")
        list = line.split(",")
    line_no += 1

para_list.close()

bind_list = " ".join(list)

time.sleep(1)

if not os.path.isfile("regi_state.txt"):
    if os.path.isfile("mean_var.txt"):
        with open("mean_var.txt", mode="r") as f:
            lnum = 0
            for line in f:
                lnum += 1
            if lnum > 0:
                time.sleep(1)
                subprocess.run("gzip -9 bin_contribution.txt", shell=True)
                exit()       

subprocess.run("rm *.out", shell=True)
subprocess.run("g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o test.out", shell=True)
subprocess.run("./test.out " + bind_list, shell=True)
subprocess.run("rm regi_state*.txt", shell=True)
subprocess.run("gzip -9 bin_contribution.txt", shell=True)

exit()