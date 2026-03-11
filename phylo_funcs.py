from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

def visualize(fasta_file):
  seqs = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
  names = [record.id for record in SeqIO.parse(fasta_file, "fasta")]

  # 塩基を数値に変換
  mapping = {"A":0,"C":1,"G":2,"T":3,"-":4}
  data = np.array([[mapping.get(base,5) for base in seq] for seq in seqs])

  # MEGA風カラー
  cmap = ListedColormap([
      "#8ED08E",  # A  緑
      "#8FA8FF",  # C  青
      "#FFD36B",  # G  オレンジ
      "#FF8A8A",  # T  赤
      "#C0C0C0",  # gap
      "#FFFFFF"
  ])

  norm = BoundaryNorm([-0.5,0.5,1.5,2.5,3.5,4.5,5.5], cmap.N)

  plt.figure(figsize=(12,2.5))
  plt.imshow(data, cmap=cmap, norm=norm, aspect="auto")

  for i in range(len(seqs)):
    for j in range(len(seqs[i])):
      plt.text(j, i, seqs[i][j], ha="center", va="center", fontsize=6)

  plt.yticks(range(len(names)), names)
  plt.xlabel("Position")
  plt.ylabel("Sequence")
  plt.show()
