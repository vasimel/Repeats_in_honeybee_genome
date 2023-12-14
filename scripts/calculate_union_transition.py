import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# path to segment bed file
# calculation of combined transitions
f = ""
df = pd.read_csv(
    f,
    names=["chr", "start", "stop", "state"],
    usecols=["state"],
    sep="\t"
).to_numpy()
d = np.zeros(
    shape=(12, 12),
    dtype=np.float32
)
for i in range(0, len(df) - 1):
    d[int(df[i][0][1:]) - 1][int(df[i + 1][0][1:]) - 1] += 1
for i in range(len(d)):
    d[i] = d[i] / d[i].sum()
x = pd.DataFrame(d,index=[i for i in range(1,13)],columns=[i for i in range(1,13)])
x.to_csv("_union_states_probabilities.csv")
sns.heatmap(
    d,
    cmap="Blues",
    annot=True,
    cbar=True,
    square=True,
    xticklabels=[f"E{i}" for i in range(1, 13)],
    yticklabels=[f"E{i}" for i in range(1, 13)]
)
plt.show()
