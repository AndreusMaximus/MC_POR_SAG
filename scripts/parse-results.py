t_state_reduction = 0
t_cpu_reduction = 0
t_mem_reduction = 0
t_false_positive = 0
t_false_negative = 0
t_correct = 0

with open('result.txt','r') as f:
    lines = f.readlines()
t_tested = 0;
for i in range(0,len(lines),2):
    t_tested += 1
    line_normal = lines[i].split(",")
    line_por    = lines[i+1].split(",")
    print(line_normal[0].split("/")[-1])

    if int(line_normal[1]) < int(line_por[1]):
        t_false_negative += 1
    if int(line_normal[1]) == int(line_por[1]):
        t_correct += 1
    if int(line_normal[1]) > int(line_por[1]):
        t_false_negative += 1

    state_reduction = (float(line_por[3])/float(line_normal[3]))*100
    cpu_reduction = (float(line_por[6])/float(line_normal[6]))*100
    mem_reduction = (float(line_por[7])/float(line_normal[7]))*100

    print("\t states reduced by: " + str(state_reduction) + "%")
    print("\t cpu reduced by   : " + str(cpu_reduction) + "%")
    print("\t mem reduced by   : " + str(mem_reduction) + "%")

    t_state_reduction += state_reduction
    t_cpu_reduction += cpu_reduction
    t_mem_reduction += mem_reduction

print("||============================results============================||")
print("|| Correct results:             " + str(t_correct))
print("|| False negatives results:     " + str(t_false_negative))
print("|| Incrorrect results:          " + str(t_false_positive))
print("|| average state reduction:     " + str(t_state_reduction/t_tested) + "%")
print("|| average cpu difference:    +-" + str(t_cpu_reduction/t_tested) + "%")
print("|| average memory difference: +-" + str(t_mem_reduction/t_tested) + "%")



