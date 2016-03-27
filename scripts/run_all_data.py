#!/usr/bin/python

import math
import os
import random
import sys

from os import listdir
from os.path import isfile, join, splitext


PROJECT_PATH = "/home/vp76/drop/cp/"
WORK_DIR = os.getcwd() + '/' #"/home/vp76/drop/4_graphs/"
EST_EXE_PATH = PROJECT_PATH + "build/indropest"
R_EXE_PATH = PROJECT_PATH + "scripts/analyze_mit.r"


def usage_exit():
    print("Usage: run_all_data.py -est | -r | -cv")
    exit(1)


def submit_r_job(hours_count, mem_gb_usage, rds_path):
    command = 'bsub -q short -W %d:00 -R "rusage[mem=%d000]" "%s %s"' % (hours_count, mem_gb_usage, R_EXE_PATH, rds_path)
    print(command)
    os.system(command)


def submit_est_job(hours_count, mem_gb_usage, bams_mask, result_name, config_path = (PROJECT_PATH + 'config.xml')):
    command = 'bsub -q short -W %d:00 -R "rusage[mem=%d000]" %s -v -m -l %s -o %s.rds -c %s %s' % (
        hours_count, mem_gb_usage, EST_EXE_PATH, result_name, WORK_DIR + result_name, config_path, bams_mask)

    print(command)
    os.system(command)


def run_all_est():
    submit_est_job(4, 45, "/groups/pklab/peterk/jimin/allon/ES/SRR1784310.tag.*.aligned.bam", 'SRR1784310.tag.full')
    submit_est_job(2, 20, "/groups/pklab/peterk/ninib/huidan/test.Dec14/hiseq/mouse/BAa1-1_*.bam", 'BAa1-1')
    submit_est_job(1, 10, "/groups/pklab/peterk/jimin/tests/K562.CK1750/SCG_29/mouse/SCG29*.aligned.bam", 'SCG29')
    submit_est_job(1, 10, "/groups/pklab/peterk/ninib/huidan/test.Dec14/mouse/MurineSample1_S2_L001.1.aligned.bam", 'MurineSample1_S2_L001.1')


def run_all_r():
    #submit_r_job(1, 7, WORK_DIR + "SCG29.rds")
    #submit_r_job(1, 7, WORK_DIR + "MurineSample1_S2_L001.1.rds")
    #submit_r_job(1, 30, WORK_DIR + "BAa1-1.rds")
    rds_names = [join(WORK_DIR, f) for f in listdir(WORK_DIR) if isfile(join(WORK_DIR, f)) and splitext(f)[1] == '.rds']
    for name in rds_names:
        submit_r_job(1, 7, name)


def run_cv_full():
    name_mask = "/groups/pklab/peterk/jimin/allon/ES/SRR1784310.tag.%d.aligned.bam"
    start_num = 1
    max_num = 27
    length = max_num - start_num

    part_size = start_num
    step_size = 1
    examples_count = 10
    jobs_count = 0
    while part_size <= max_num:
        for i in range(examples_count):
            jobs_count += 1
            bam_str = ''.join([name_mask % rn + ' ' for rn in random.sample(range(1, max_num + 1), part_size)])
            submit_est_job(1 + int(part_size / 9), 3 + int((part_size - 1) * 1.5) , bam_str, 'SRR1784310_cv_%d_%d' % (part_size, i))

        part_size += step_size

    print(jobs_count)


def run_cv_3():
    name_mask = "/groups/pklab/peterk/jimin/allon/ES/SRR1784310.tag.%d.aligned.bam"
    max_num = 27

    for part_size in range(1, max_num + 1):
        bam_str = ''.join([name_mask % random.randint(1, max_num) + ' ' for _ in range(1, int(part_size) + 1)])
        submit_est_job(1 + int(part_size / 9), 3 + int((part_size - 1) * 1.5) , bam_str, 'SRR1784310_cv_%d_1' % (part_size))


def run_cv():
    name_mask = "/groups/pklab/peterk/jimin/allon/ES/SRR1784310.tag.%d.aligned.bam"
    start_num = 2
    max_num = 27
    length = max_num - start_num

    part_size = start_num
    step_size = 1.5
    examples_count = 10
    examples_step_size = examples_count * math.log(step_size) / math.log(length)
    jobs_count = 0
    while part_size <= max_num:
        print(str(int(examples_count) - 1) + ':')
        for i in range(1, int(examples_count)):
            jobs_count += 1
            bam_str = ''.join([name_mask % random.randint(1, max_num) + ' ' for _ in range(1, int(part_size) + 1)])
            submit_est_job(1 + int(part_size / 9), 3 + int(part_size - 1) * 2, bam_str, 'SRR1784310_cv_%d_%d' % (part_size, i))

        examples_count -= examples_step_size

        part_size *= step_size

    print(jobs_count)


def run_pred_stat():
    name_mask = "/groups/pklab/peterk/jimin/allon/ES/SRR1784310.tag.%d.aligned.bam"
    for ex in range(100):
        bam_str = ''.join([name_mask % id + ' ' for id in random.sample(range(1,28), 5)])
        submit_est_job(1, 11, bam_str, 'SRR1784310_cv_%d_0' % ex)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        usage_exit()

    if sys.argv[1] == '-est':
        run_all_est()
    elif sys.argv[1] == '-r':
        run_all_r()
    elif sys.argv[1] == '-cv':
        run_cv()
    elif sys.argv[1] == '-cv2':
        run_cv_full()
    elif sys.argv[1] == '-cv3':
        run_cv_3()
    elif sys.argv[1] == '-ps':
        run_pred_stat()
    else:
        usage_exit()
