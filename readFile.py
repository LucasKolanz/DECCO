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

    for line in content.split("step: ")[1:-1]:
        elements = line.split("\n")
        steps.append(int(elements[0]))
        totalTorque.append(float(elements[1][len("totalTorqueA: 0,0,"):]))
        w.append(float(elements[4][len("w[1]: 0,0,"):]))
        nA.append(float(elements[5][len("nA: "):].split(',')[0]))
        power.append(float(elements[6][len("power: "):]))

    print(f"Time Integrated Power: {np.trapz(power)}")
    start = 0
    end = 10000

    fig,ax = plt.subplots()

    maxTorque = max([abs(i) for i in totalTorque])
    maxw = max([abs(i) for i in w])
    maxnA = max([abs(i) for i in nA])

    maxmax = max(maxTorque,maxw,maxnA)
    # maxmax = 1

    totalTorque = [i*(maxmax/maxTorque) for i in totalTorque]
    w = [i*(maxmax/maxw) for i in w]
    nA = [i*(maxmax/maxnA) for i in nA]

    print(f"min torque in range {start} to {end}: {min(totalTorque[start:end])} at step {start + totalTorque[start:end].index(min(totalTorque[start:end]))}")
    print(f"min w in range {start} to {end}: {min(w[start:end])} at step {start + w[start:end].index(min(w[start:end]))}")
    print(f"min nA in range {start} to {end}: {min(nA[start:end])} at step {start + nA[start:end].index(min(nA[start:end]))}")

    ax.plot(steps[start:end],totalTorque[start:end],label="totalTorque",marker=".")
    ax.plot(steps[start:end],w[start:end],label=f"w",marker="*")
    ax.plot(steps[start:end],nA[start:end],label=f"nA")

    ax.legend()

    plt.show() 

    # for data in [totalTorque,w,nA]:


    #     fig,ax = plt.subplots()


    #     ax.plot(steps[:end],data[:end])
    #     # ax.plot(steps[:end],w[:end],label=f"w")
    #     # ax.plot(steps[:end],nA[:end],label=f"nA")

    #     plt.show() 
