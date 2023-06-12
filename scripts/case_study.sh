#!/bin/bash
rm case_study.csv
echo "j0"
echo "SAG"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j0-c0.csv

echo "SAG-POR"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j0-c0.csv --por=priority

echo "SAG-POR-interfering"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j0-c0.csv --por=priority --interfering=all

echo "SAG-POR-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j0-c0.csv --por=priority --limit_fail=yes

echo "SAG-POR-interfering-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j0-c0.csv --por=priority --interfering=all --limit_fail=yes

echo "j0"
echo "SAG"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j5-c0.csv

echo "SAG-POR"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j5-c0.csv --por=priority

echo "SAG-POR-interfering"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j5-c0.csv --por=priority --interfering=all

echo "SAG-POR-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j5-c0.csv --por=priority --limit_fail=yes

echo "SAG-POR-interfering-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j5-c0.csv --por=priority --interfering=all --limit_fail=yes

echo "j10"
echo "SAG"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j10-c0.csv

echo "SAG-POR"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j10-c0.csv --por=priority

echo "SAG-POR-interfering"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j10-c0.csv --por=priority --interfering=all

echo "SAG-POR-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j10-c0.csv --por=priority --limit_fail=yes

echo "SAG-POR-interfering-limit-fails"
./../build/nptest -m 4 ../../Re\ Experiment\ Results/mer-app-frame-job-j10-c0.csv --por=priority --interfering=all --limit_fail=yes
