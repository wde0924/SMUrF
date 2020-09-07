import os
import time
import glob
import datetime
import subprocess

#cmd_template = ["Rscript", "--vanilla", "main_script_Reco.R"]
cmd_template = ["Rscript", "--vanilla", "main_script_NEE.R"]

t = 60 * 5  # sec
pid = None
idx = 0
jobs = []

#'westernCONUS',    'easternCONUS', 'westernEurope', 'easternChina', 
#'easternAustralia' 'easternAsia'   'southAmerica'

reg_indx = [1, 2, 3]
#yrs = ['2010', '2011', '2012', '2013', '2014']
yrs = ['2013']
mons = range(7, 12)

for r in reg_indx:
    for y in yrs:
        for m in mons: jobs.append([r, y, '{:02d}'.format(m)])

print("Generating jobs[{}]: {}".format(len(jobs), jobs))
#exit()

start_time = None
end_time = None

def submit():
    global idx
    global pid
    global start_time

    cmd = cmd_template + list(map(str, jobs[idx]))
    job = subprocess.run(cmd, stdout=subprocess.PIPE)
    job_stdout = str(job.stdout)

    # search for PID after "submitted batch job" from running the main script
    pid_line_idx = job_stdout.find("Submitted batch job")
    if pid_line_idx != -1:
        pid = job_stdout[pid_line_idx:].split(" ")[-1][:-3]
        start_time = datetime.datetime.now()
        print("Caught {} at {} ".format(pid, start_time))
    else:
        print("Submit job failed!")
        print(" ".join(cmd))
    idx = idx + 1

submit()

while True:
    if pid == None:
        print("No job is currently running.")
        submit()
        continue

    current_jobs = str(subprocess.run(["squeue", "-A", "lin-kp"], stdout=subprocess.PIPE).stdout)

    if current_jobs.find(pid) == -1:
        end_time = datetime.datetime.now()
        print("{} had finished at {}".format(pid, end_time))
        print("{} time diff: {}".format(pid, end_time - start_time))
        if idx >= len(jobs):
            print("all jobs had finished")
            break
        else:
            submit()
    else:
        print("{} is still running".format(pid))

    time.sleep(t)

