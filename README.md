# ML_hw
HW #1
We are going to use pre-processed RNA-seq data from mouse brain cells with 2 genotypes - wild type (WT) and ADAR1-KO (aka knockout/inactive/disrupted ADAR1 gene).
```
import pandas as pd
import numpy as np

# Load tab separated HW-1.TPM.tsv.gz file, and create a data frame 
# hint: read_csv with sep
df = pd.read_csv('/content/HW-1.TPM.tsv.gz',sep='\t')

# Show the first 5 rows
df.head()
# Before continuing, we also need to check that the numeric columns have been parsed correctly and that we don't have any missing values.
# Check that we don't have missing values
# hint: isna with sum
totalna = df.isna().sum()
print("Total NA cells: ", totalna)

# Assert that its actually 0 and provide some message if it's not
assert df.any != 'na', 'Total NA cells !=0'
# The values ​​of the normalized expression in this dataset cannot be less than 0. Let's check if it's true:
# Select numerical columns
# hint: select_dtypes
numcols = df.select_dtypes(include=np.number).columns
# Assert that the total number of values < 0 is 0
f=(df[numcols]<0).sum()
assert f.sum() == 0, \
  "All RNA abundance estimates must be > 0"
  # Calculate non-zero quantiles for each expression column
qthr = 0.01
quantiles = []
for col in df[numcols]:
  series = df[col]
  # Select non zero values
  series = series[series>0]
  # Calculate quantiles
  q = series.quantile(qthr)
  quantiles.append(q)

# Print results
for col, q in zip(numcols, quantiles):
  print(df[col].name,'<-',q)

  # assert that quantile is not zero
  assert q != 0, \
  f*quantile (df[col].name)
  # As a threshold we will use a min quantile
threshold = min(quantiles)

print("Genes before", df.shape[0])

# We will drop all columns, 
# where the expression is below the threshold in all samples
mask = (df[numcols] >= threshold).any(axis=1)
df = df[mask]

print("Genes after", df.shape[0])
# Let's check the genes with the total maximum expression in all samples:
overallexpr = df[numcols].sum(axis=1)

# New trick:
## sort by index
argsort = overallexpr.argsort()
## select top 25 elements
index = argsort[-15:]

# Print selected genes
# hint: use iloc for indexing using row ids
df.iloc[index]
# It is a well-known observation that PCA decomposition of expression profiles should group samples according to their origin. In our case, we expect to observe two clusters - KO and WT cells.
# The closer the samples are to each other on the PCA plot, the better. However, we are fine as long as the samples can be separated by a straight line.

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Select expression columns 
data = df[numcols]
# Transpose the matrix to treat genes as features
data = data.transpose()
# Transform to zero mean and unit variance
data = StandardScaler().fit_transform(data)

# first thing is to calculate the PCA decomposition
# hint: don't forget the random state & request 2 components to make a 2D plot
pca = PCA(random_state=32, n_components=2).fit_transform(data)

samples, coords = pca.shape
print(f"Samples: {samples}")
print(f"Components: {coords}")
import matplotlib.pyplot as plt
# colors in HEX format
palette = {"WT": "#1F77B4", "ADAR1-KO": "#D62728"}

# Create a basic draft
fig = plt.figure()
ax = fig.gca()

for (x, y), col in zip(pca, numcols):
  group = "WT" if "WT" in col else "ADAR1-KO"
  color = palette[group]
  ax.scatter(x, y, color=color)

fig.show()
palette = {"WT": "#1F77B4", "ADAR1-KO": "#D62728"}

# Create a basic draft
fig = plt.figure()
ax = fig.gca()

for (x, y), col in zip(pca, numcols):
  group = "WT" if "WT" in col else "ADAR1-KO"
  color = palette[group]
  ax.scatter(x, y, color=color)

  plt.grid()
  ax.set_xlabel('PCA 1', fontstyle='italic')
  ax.set_ylabel('PCA 2', fontstyle='italic')
  plt.title("PCA decomposition", fontweight='bold', loc='left', fontsize=25)
  plt.text(-130,120, "ADAR1 KO", color = "#D62728", fontweight="bold", fontsize=15)
  plt.text(110,120, "WT", color = "#1F7784", fontweight="bold", fontsize=15)

#fig.show()
# We can make a more formal version using the distance matrix (calculated over expression values).
import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

# Calculate the spearman correlation
p_data=df[numcols]
corr = p_data.corr(method='spearman')
corr = 1 - corr


fig, ax = plt.subplots() # same as figure + .gca()

# cmap - mapping between values and colors
cmap = cmap=plt.get_cmap('hot')
# imshow = image show (yes, we treat the matrix as a picture)
im = ax.imshow(corr, cmap=cmap)

loc, labels = np.arange(len(numcols)), list(numcols)

# Create ticks
ax.set_xticks(loc)
ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")

ax.set_yticks(loc)
ax.set_yticklabels(labels, fontsize=13)

# label cells
for x in loc:
  for y in loc:
    value = corr.iloc[x, y]
    ax.text(x-0.25,y,round(value,3),fontsize=12)

# Disable spines
for s in ax.spines.values():
  s.set_visible(False)

# Add minor ticks
ax.set_xticks(loc - 0.5, minor=True)
ax.set_yticks(loc - 0.5, minor=True)

# And create a grid based on them
ax.grid(which='minor',color='white',linewidth=4)
ax.tick_params(which="minor", bottom=False, left=False)

ax.set_title("Expr.correlation", loc='left',fontweight="bold")

# Colorbar to visualize the cmap
fig.colorbar(im, ax=ax)

fig.show()
# Let's have a look at binned distribution of genes expression:
fig, ax = plt.subplots()

# Let's ignore the highly expressed genes - 
# it makes no sense to bin them, since there are only a few of them.
upthr = df[numcols].quantile(0.98).max()

# Create 200 bins ranging from 0 to upthr
# hint: use numpy
bins = np.arange(0,upthr,upthr/200)


# Calculate center of each bin
# hint: slice with step & sum & divide by 2
X = [(bins[i-1]+bins[i])/2 for i in range(1,200)]

# Plot each genotype
for genotype in "WT", "ADAR1-KO":

  # There are 3 samples of each genotype
  counts = []
  for col in numcols:
    if genotype not in col:
      continue
    # count number of genes in each bin
    c, _ = np.histogram(df[col], bins)
    counts.append(c)
  
  # Turn list of 1D arrays into dense 2D array
  counts = np.asarray(counts)
  # Calculate stats for each bin
  meanv, minv, maxv = counts.mean(axis=0), counts.min(axis=0), counts.max(axis=0)

  # Get the color from the palette
  color = 'green'

  # Plot the mean trend
  ax.plot(X,meanv)
  # Shade betweenmin and max
  ax.fill_between(X,minv,maxv, color=color, alpha=0.25)

ax.set(yscale='log')
# Now your task is to use all your skills and embellish this plot just like we did with the PCA figure. For example, be sure to include the legend, title, grid, and axis labels.
from matplotlib.pyplot import title,xlabel,ylabel
fig, ax = plt.subplots()
upthr = df[numcols].quantile(0.98).max()

bins = np.arange(0,upthr,upthr/200)

X = [(bins[i-1]+bins[i])/2 for i in range(1,200)]

for genotype, color_line in zip(["WT", "ADAR1-KO"],['green','red']):
  counts = []
  for col in numcols:
    if genotype not in col:
      continue
    c, _ = np.histogram(df[col], bins)
    counts.append(c)
 
  counts = np.asarray(counts)
  meanv, minv, maxv = counts.mean(axis=0), counts.min(axis=0), counts.max(axis=0)

  ax.plot(X,meanv,color=color_line)
  # Shade betweenmin and max
  ax.fill_between(X,minv,maxv, color='dark' + color_line, alpha=0.25)
  
  ax.set(yscale='log', ylabel='expressed genes', xlabel='expression')
  ax.legend(title='expr.distribution')
  ax.legend(['WT','WT min-max','ADAR1-KO','ADAR1-KO min-max'])
  ax.grid()

ax.set(yscale='log')
# Now let's visualize the expression of top N genes. First - let's get the data:
# As before, we will select genes based on their overall expression
sumexpr = data.sum(axis=0).argsort()

topn = 250
# select indices of topn genes with max expression
ind = sumexpr[-topn:]

# get the expression values
expr = data[:,ind]
expr

import seaborn as sns

# Plot the data
cluster = expr.transpose()

# Disable x axis labels
# hint: cluster.ax_heatmap is a simple matplotlib axis instance
map=sns.clustermap(cluster, yticklabels=False)
map.fig.suptitle('Expression of top 250 genes', fontsize=15,fontweight='bold')
```
HW #2
```
#скачиваем нужные данные
# 1 ChIP-seq replica
!wget -O A549-1.bigBed "https://www.encodeproject.org/files/ENCFF147PDF/@@download/ENCFF147PDF.bigBed"
!wget -O A549-2.bigBed "https://www.encodeproject.org/files/ENCFF479CBZ/@@download/ENCFF479CBZ.bigBed"
# 3 ATAC-seq replicas
!wget -O Chip-seq1.bigBed "https://www.encodeproject.org/files/ENCFF297TBO/@@download/ENCFF297TBO.bigBed"
!wget -O Chip-seq2.bigBed "https://www.encodeproject.org/files/ENCFF120BYY/@@download/ENCFF120BYY.bigBed"
!wget -O Chip-seq3.bigBed "https://www.encodeproject.org/files/ENCFF352AOI/@@download/ENCFF352AOI.bigBed"

!wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
!chmod a+x bigBedToBed
# конверируем файлы в нужный формат
for file in "A549-1","A549-2", "Chip-seq1", "Chip-seq2", "Chip-seq3":
  !./bigBedToBed "{file}.bigBed" "{file}.bed"

!apt install -y bedtools
!pip3 install pybedtools

from pybedtools import BedTool
#чтобы дальше было удобнее работать
ATACseq1 = BedTool("A549-1.bed").sort()
ATACseq2 = BedTool("A549-2.bed").sort()
chip1 = BedTool("Chip-seq1.bed").sort()
chip2 = BedTool("Chip-seq2.bed").sort()
chip3 = BedTool("Chip-seq3.bed").sort()

# уникальные пики для ATACseq пересечений
ATAC = ATACseq1.intersect(ATACseq2).sort()

ATACseq1_only = ATACseq1.subtract(ATACseq2).sort()
ATACseq2_only = ATACseq2.subtract(ATACseq1).sort()

atac_not_replicated = ATACseq1_only.cat(ATACseq2_only).sort()

# проверка
assert ATAC.intersect(atac_not_replicated).total_coverage() == 0

#удалим области из данных Chip-seq
chip1 = chip1.subtract(atac_not_replicated).sort()
chip2 = chip2.subtract(atac_not_replicated).sort()
chip3 = chip3.subtract(atac_not_replicated).sort()
# -wa means "keep a whole A549 peak in case of an overlap with ATAC-seq"
fg1 = chip1.intersect(ATAC, wa=True, u=True).sort()
fg2 = chip2.intersect(ATAC, wa=True, u=True).sort()
fg3 = chip3.intersect(ATAC, wa=True, u=True).sort()
# -A means "remove entire ATAC peaks overlapping ATAC-seq peak"
bg1 = ATAC.subtract(chip1, A=True).sort()
bg2 = ATAC.subtract(chip2, A=True).sort()
bg3 = ATAC.subtract(chip3, A=True).sort()

# проверка
assert fg1.intersect(bg1).total_coverage() == 0
assert fg2.intersect(bg2).total_coverage() == 0
assert fg3.intersect(bg3).total_coverage() == 0

TF_list = (chip1,chip2,chip3,ATAC)

TF_classes = ["chip1","chip2","chip3","ATAC"]
seq = [chip1,chip2,chip3,ATAC]

for i,tf_i in enumerate(seq):
   tbl_row = []
   for tf_j in seq:
     tbl_row.append(
         round(len(tf_i.intersect(tf_j, wa = True, u = True).sort())/
                   len(tf_i), 2)
     )
   print(TF_classes[i], end = '\t')
   print(*tbl_row, sep = '\t')
   
print('\n', end = '\t')
print(*TF_classes, sep = '\t')

OHE = pd.get_dummies(TF_classes)
OHE

import seaborn as sns
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(20, 10))

for title, regions, bw, ax in \
    ("SMC3 & ATAC-seq", fg1, 25, axes[0,0]), \
    ("RAD21 & ATAC-seq", fg2, 25, axes[0,1]), \
    ("EHF & ATAC-seq", fg3, 25, axes[1,0]), \
    ("ATAC-seq", bg1, 100, axes[1,1]):
  # гистограмма
  sns.histplot([x.length for x in regions], binwidth=bw, kde=True, ax=ax)
  ax.set(title=title, xlabel="Peak length", ylabel="Density")

# получение ДНК последовательностей: 
!gsutil -m cp \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
  .

fasta = "Homo_sapiens_assembly38.fasta"

# извлечем таргетные последовательности
# seqfn - path to final FASTA file
fgseq1 = fg1.sequence(fi=fasta).seqfn
fgseq2 = fg2.sequence(fi=fasta).seqfn
fgseq3 = fg3.sequence(fi=fasta).seqfn
bgseq = bg1.sequence(fi=fasta).seqfn

!pip3 install biopython
from Bio import SeqIO

# парсинг
fgseq1 = [str(x.seq) for x in SeqIO.parse(fgseq1, format='fasta')]
fgseq2 = [str(x.seq) for x in SeqIO.parse(fgseq2, format='fasta')]
fgseq3 = [str(x.seq) for x in SeqIO.parse(fgseq3, format='fasta')]

bgseq = [str(x.seq) for x in SeqIO.parse(bgseq, format='fasta')]

# проверка
print("SMC3", fgseq1[0])
print("RAD21", fgseq2[0])
print("EHF", fgseq3[0])
import numpy as np

print(f"Before \t fg: {len(fgseq1)}; bg: {len(bgseq)}")
print(f"Before \t fg: {len(fgseq2)}; bg: {len(bgseq)}")
print(f"Before \t fg: {len(fgseq3)}; bg: {len(bgseq)}")
#так как это учебное задание, можно сократить кол-во данных для более быстрого получения результата

np.random.seed(123)
fgseq1 = np.random.choice(fgseq1, size=2_000, replace=False)
fgseq2 = np.random.choice(fgseq2, size=2_000, replace=False)
fgseq3 = np.random.choice(fgseq3, size=2_000, replace=False)
bgseq = np.random.choice(bgseq, size=8_000, replace=False)

print(f"After \t fg: {len(fgseq1)}; bg: {len(bgseq)}")
print(f"After \t fg: {len(fgseq2)}; bg: {len(bgseq)}")
print(f"After \t fg: {len(fgseq3)}; bg: {len(bgseq)}")

#сделаем фичи
# подсчитаем k-меры
from collections import defaultdict

def calculate_kmers(seq: str, klen: int):
  assert len(seq) >= klen and klen >= 1, seq
  total_kmers = len(seq) - klen + 1

  counts = defaultdict(int)
  for ind in range(total_kmers):
    window = seq[ind:ind+klen]
    counts[window] += 1

  # исключим non-ATGC k-mers
  counts = {
      k: v for k, v in counts.items() if {"A", "C", "G", "T"}.issuperset(set(k))
  }

  # подсчитаем частоты
  total_kmers = sum(counts.values())
  frequencies = {k: v / total_kmers for k, v in counts.items()}
  return frequencies

from tqdm import tqdm

KMERS = 1, 2, 3, 4, 5

# будущий DataFrame
df = []
for cls, sequences in (0, bgseq), (1, fgseq1):
  # tqdm draws progress bar while iterating over collection
  for seq in tqdm(sequences):
    record = {}
    for klen in KMERS:
      record.update(calculate_kmers(seq, klen))
    record['Class'] = cls
    df.append(record)

import pandas as pd

df = pd.DataFrame(df).fillna(0)
df.head()

4**5 + 4**4 + 4**3 + 4**2 + 4**1 + 1

# проверка для labels
df['Class'].value_counts()

# проверка для dtypes
df.dtypes.value_counts()

df.describe()

#теперь поделим данные на трейн и тест
Y = df.pop('Class').values
features = df.columns.values
X = df.values

from sklearn.model_selection import train_test_split

Xtrain, Xtest, Ytrain, Ytest = train_test_split(
    X, Y, test_size=0.3, random_state = 123
)

print("Train:")
print(f"\tX: {Xtrain.shape}; Y: {Ytrain.shape}")
print("Test:")
print(f"\tX: {Xtest.shape}; Y: {Ytest.shape}")

# Logistic regression 

from sklearn.linear_model import LogisticRegression
from sklearn.multiclass import OneVsRestClassifier

# Определяем набор данных

#так как у нас несколько классов, то нам нужно использовать One-vs-rest стратегию
model = LogisticRegression()
ovr = OneVsRestClassifier(model)
# фитим модель нашими данными
ovr.fit(Xtrain, Ytrain)
# и делаем предсказание
Ypred = ovr.predict(Xtrain)
Ytrue = Ytrain

# смотрим что предсказалось на 50 образцах
print(f"Predicted:\n {Ypred[:20]}")
print(f"Expected:\n {Ytrue[:20]}")

# как можно видеть, чаще всего модель говорит, что у нас объект 0 класса, при этом объекты других классов не предсказываются
# посчитаем accuracy

from sklearn.metrics import accuracy_score
accuracy_score(Ytrue, Ypred)

# как можно видеть, accuracy 0,798. То есть наша модель предсказывает около 80% данных.
# посмотрим другие метрики

# для начала посчитаем, как точно модель предсказывает положительные классы
from sklearn.metrics import precision_score
ps = precision_score(Ytrue, Ypred, average = 'micro')
print(f"Precision score: {ps: .3f}")
# не так плохо. Посмотрим для каждого класса:
ps = precision_score(Ytrue, Ypred, average = None)
print(ps)
# как можно видеть модель сравнительно неплохо предсказывает 0 класс
# в принципе уже понятно, что линейная регрессия не годится, но посчитаем для практики другие метрики

# сколько положительных были предсказаны моделью как положительные
from sklearn.metrics import recall_score
recall_score(Ytrue, Ypred, average = 'micro')
#больше половины

# оценим эффективность нашего классификатора в целом
# F1 будет более показателен, так как оно включает в себя точность и полноту наших предсказаний(precision и recall)
from sklearn.metrics import f1_score
f1_score(Ytrue, Ypred, average='micro')
# В результате удается предсказать 80% наших данных - не совсем провал

from sklearn.metrics import classification_report
rep = classification_report(Ytrue, Ypred)
print(rep)

from sklearn.metrics import roc_auc_score
roc_auc_score(Ytrue, Ypred, average=None)

#Decision Tree

from sklearn.tree import DecisionTreeClassifier

tree = DecisionTreeClassifier(random_state=17, max_depth=10).fit(Xtrain, Ytrain)
Y_pred = tree.predict(Xtrain)

rep = classification_report(Ytrue, Ypred)
print(rep)

from sklearn.model_selection import GridSearchCV

grid = {
    'max_depth': [3, 6, 9, 12],
    'class_weight': [None, 'balanced'],
    'min_samples_leaf': [50, 100, 200]
}

grid = GridSearchCV(
    DecisionTreeClassifier(), grid, cv=3, scoring='f1_micro'
).fit(Xtrain, Ytrain)

Y_pred = tree.predict(Xtrain)

rep = classification_report(Ytrue, Ypred)
print(rep)

#random forest
from sklearn.ensemble import RandomForestClassifier

kwargs = dict(
    max_depth=6, min_samples_leaf=50, n_jobs=-1, 
    class_weight='balanced_subsample'
)

forest = RandomForestClassifier(n_estimators=50, **kwargs).fit(Xtrain, Ytrain)
print("Random forest (N=50):", tree, Xtest, Ytest)

#SVM

from sklearn.svm import SVC

svm = SVC(probability=True)
svm.fit(Xtrain, Ytrain)

print("Naive SVM:", svm, Xtest, Ytest)
