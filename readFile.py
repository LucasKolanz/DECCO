import matplotlib.pyplot as plt

import numpy as np



def read_entire_file(filepath):
    try:
        with open(filepath, 'r') as file:
            content = file.read()
            return content
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
    except IOError as e:
        print(f"I/O error: {e}")

if __name__ == "__main__":
    
    content = read_entire_file("/mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_branch/SpaceLab_data/jobs/JKRTest/test.txt")

    steps = []
    totalTorque = []
    w = []
    nA = []
    power = []
    vel = []
    rollingDisp = []
    slidingDisp = []    

    for line in content.split("step: ")[1:-1]:
        elements = line.split("\n")
        steps.append(int(elements[0]))
        totalTorque.append(float(elements[1][len("totalTorqueA: 0,0,"):]))
        w.append(float(elements[4][len("w[1]: 0,0,"):]))
        nA.append(float(elements[5][len("nA: "):].split(',')[1]))
        power.append(float(elements[6][len("power: "):]))

        vel.append(float(elements[7][len("VelNorm: "):]))
        rolldisp = elements[8][len("rollingDisp: "):].split(',')
        rollingDisp.append(np.linalg.norm([float(rolldisp[0]),float(rolldisp[1]),float(rolldisp[2])]))
        sliddisp = elements[9][len("slidingDisp: "):].split(',')
        # slidingDisp.append(np.linalg.norm([float(sliddisp[0]),float(sliddisp[1]),float(sliddisp[2])]))
        slidingDisp.append(float(sliddisp[0]))

    # print(f"Time Integrated Power: {np.trapz(power)}")
    start = 1
    end = -1
    step = 1

    fig,ax = plt.subplots()

    # maxTorque = max([abs(i) for i in totalTorque])
    # maxw = max([abs(i) for i in w])
    # # maxw=0
    # maxnA = max([abs(i) for i in nA])
    # maxmax = max(maxTorque,maxw,maxnA)

    # # maxmax = 1

    # totalTorque = [i*(maxmax/maxTorque) for i in totalTorque]
    # w = [i*(maxmax/maxw) for i in w]
    # # nA = [i*(maxmax/maxnA) for i in nA]

    print(f"min torque in range {start} to {end}: {min(totalTorque[start:end])} at step {start + totalTorque[start:end].index(min(totalTorque[start:end]))}")
    print(f"corresponding w in range {start} to {end}: {w[totalTorque[start:end].index(min(totalTorque[start:end]))]} at step {start + totalTorque[start:end].index(min(totalTorque[start:end]))}")
    print(f"min nA in range {start} to {end}: {min(nA[start:end])} at step {start + nA[start:end].index(min(nA[start:end]))}")

    # ax.plot(steps[start:end],totalTorque[start:end],label="totalTorque",marker=".")
    ax.plot(steps[start:end],w[start:end],label=f"w",marker="*")
    ax.plot(steps[start:end],nA[start:end],label=f"nA")
    # ax.plot(steps[start:end],vel[start:end],label=f"vel")
    # ax.plot(list(range(len(vel[start:end]))),vel[start:end],label=f"vel")
    # ax.plot(steps[start:end],rollingDisp[start:end],label=f"rollingDisp")
    ax.plot(steps[start:end],slidingDisp[start:end],label=f"slidingDisp")
    # ax.axhline(y=0)

    ax.legend()

    plt.show() 

    # fig,ax = plt.subplots()
    # ax.plot(steps[start:end:step],power[start:end:step],label="power",marker=".")
    # plt.show() 
    # for data in [totalTorque,w,nA]:


    #     fig,ax = plt.subplots()


    #     ax.plot(steps[:end],data[:end])
    #     # ax.plot(steps[:end],w[:end],label=f"w")
    #     # ax.plot(steps[:end],nA[:end],label=f"nA")

    #     plt.show() 
