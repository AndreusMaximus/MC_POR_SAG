import sys

t_state_reduction = 0
t_cpu_reduction = 0
t_mem_reduction = 0
t_false_positive = 0
t_false_negative = 0
t_correct = 0

t_s_reductions = 0
t_f_reductions = 0

t_sched_normal = 0
t_sched_por = 0

t_timeout_normal = 0
t_timeout_por = 0

with open('result.txt','r') as f:
    lines = f.readlines()
t_tested = 0
for i in range(0,len(lines),2):
    t_tested += 1
    line_normal = lines[i].split(",")
    line_por    = lines[i+1].split(",")
    #print(line_normal[0].split("/")[-1])
    t_timeout_normal += int(line_normal[8])
    t_timeout_por += int(line_por[8])

    t_sched_normal += int(line_normal[1])
    t_sched_por += int(line_por[1])

    state_reduction = (float(line_por[3])/float(line_normal[3]))*100
    cpu_reduction = (float(line_por[6])/float(line_normal[6]))*100
    mem_reduction = (float(line_por[7])/float(line_normal[7]))*100

    if int(line_normal[1]) < int(line_por[1]):
        print(line_normal[0])
        t_false_positive += 1
    if int(line_normal[1]) == int(line_por[1]):
        t_correct += 1
        t_state_reduction += state_reduction
        t_cpu_reduction += cpu_reduction
        t_mem_reduction += mem_reduction
    if int(line_normal[1]) > int(line_por[1]):
        t_false_negative += 1


    #print("\t states reduced by: " + str(state_reduction) + "%")
    #print("\t cpu reduced by   : " + str(cpu_reduction) + "%")
    #print("\t mem reduced by   : " + str(mem_reduction) + "%")

    t_s_reductions += int(line_por[10])
    t_f_reductions += int(line_por[11])

print("||============================results============================||")
print("|| cores:                       " + sys.argv[1])
print("|| utilization:                 " + sys.argv[2])
print("|| Correct results:             " + str(t_correct))
print("|| False negatives results:     " + str(t_false_negative))
print("|| Incorrect results:           " + str(t_false_positive))
print("|| average state reduction:     " + str(t_state_reduction/t_correct) + "%")
print("|| average cpu difference:    +-" + str(t_cpu_reduction/t_correct) + "%")
print("|| average memory difference: +-" + str(t_mem_reduction/t_correct) + "%")
print("|| total successfull reductions " + str(t_s_reductions))
print("|| total failed reductions      " + str(t_f_reductions))
print("|| timeouts normal - por:       " + str(t_timeout_normal) + " | " + str(t_timeout_por))
print("|| schedulable normal - por:    " + str((t_sched_normal/(t_tested))*100) + "% | " + str((t_sched_por/(t_tested))*100) + " %")



with open('out.txt','a+') as f:
    f.write("||============================results============================||" + "\n")
    f.write("|| cores:                       " + sys.argv[1] + "\n")
    f.write("|| utilization:                 " + sys.argv[2] + "\n")
    f.write("|| Correct results:             " + str(t_correct) + "\n")
    f.write("|| False negatives results:     " + str(t_false_negative) + "\n")
    f.write("|| Incorrect results:           " + str(t_false_positive) + "\n")
    f.write("|| average state reduction:     " + str(t_state_reduction/t_correct) + "%" + "\n")
    f.write("|| average cpu difference:    +-" + str(t_cpu_reduction/t_correct) + "%" + "\n")
    f.write("|| average memory difference: +-" + str(t_mem_reduction/t_correct) + "%" + "\n")
    f.write("|| total successfull reductions " + str(t_s_reductions) + "\n")
    f.write("|| total failed reductions      " + str(t_f_reductions) + "\n")
    f.write("|| timeouts normal - por:       " + str(t_timeout_normal) + " | " + str(t_timeout_por) + "\n")
    f.write("|| schedulable normal - por:    " + str((t_sched_normal/(t_tested))*100) + "% | " + str((t_sched_por/(t_tested))*100) + " %" + "\n")

