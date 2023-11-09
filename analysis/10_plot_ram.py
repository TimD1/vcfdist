import matplotlib.pyplot as plt
plt.figure(figsize=(15,5))

run = 7

# plot RAM usage
ram_file = open(f"../src/ram{run}.txt", 'r')
step = 10
ram_usage = []
time = 0
times = []
for line in ram_file:
    if line[:4] == "Mem:":
        fields = line.strip().split()
        ram = int(fields[2])
        ram_usage.append(ram / (1024*1024))
        times.append(time)
        time += step
    else:
        continue

# plot vcfdist pipeline
stages = ["read", "cluster query", "realign query", "recluster query", 
    "cluster truth", "realign truth", "recluster truth", "supercluster", 
    "precision-recall", "distance", "phase", "write"]
# stage_runtimes = [24, 311, 44, 282, 331, 44, 262, 0, 281, 0, 5] # run 0
# stage_runtimes = [24, 268, 44, 236, 308, 44, 236, 0, 281, 0, 5] # run 1
# stage_runtimes = [29, 212, 164, 236, 202, 164, 236, 0, 1718, 203, 0, 197] # run 6
stage_runtimes = [205, 518, 1159, 419, 518, 1159, 419, 0, 7939, 655, 0, 200] # run 7
colors = 'rcygbmk'
time = 0
i = 0
for stage, runtime in zip(stages, stage_runtimes):
    plt.axvspan(time, time+runtime, facecolor=colors[i%len(colors)], alpha=0.5)
    plt.text(time, 0.5 + i*2, stage)
    time += runtime
    i += 1

plt.plot(times, ram_usage, 'k.')
plt.xlabel("Time (seconds)")
plt.ylabel("RAM Usage (GiB)")
plt.ylim(0, max(ram_usage)*1.1)
plt.savefig(f"ram{run}.png", dpi=200)
