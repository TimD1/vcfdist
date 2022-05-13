import matplotlib.pyplot as plt
import numpy as np

gap_size = [0, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200]
edit_dist = [12_274_744, 12_240_524, 12_238_398, 12_235_566, 12_231_452, 12_228_310, 12_226_298, 12_225_198, 12_224_610, 12_224_272, 12_224_140, 12_224_072, 12_224_000, 12_223_803]

plt.ylabel('Edit Distance')
plt.xlabel('Gap Size for Independence')
plt.plot(gap_size, edit_dist)
plt.tight_layout()
plt.savefig('img/gap_edits.png')
