# Phylogenetics Tutorial

1. Copy phylo data

```
cp -r /analyses/vdworkshop/phylogenetics/ ./
```

2. activate mafft module 


```
module miniconda && conda activate mafft
```

3. Check mafft manual:

More details here: https://mafft.cbrc.jp/alignment/server/index.html


```
mafft --help
```


4. Examples of running:

```
#default (fastest)

mafft input.fst > input.mafft.fst

#special options (slower, check link in step 3)

mafft --localpair input.fst > input.linsi.fst

mafft --genafpair input.fst > input.einsi.fst

mafft --globalpair input.fst > input.ginsi.fst


#reordering based on sequence similarity

mafft --reorder input.fst > input.mafft.fst


```



5. Use `scp` to download the alignment from the server to your computer
and view the alignment in Aliview


6. If it looks good then go ahead with using the alignment to make a phylogeny! 


7. Check the `iqtree` manual first

```
iqtree -h
```