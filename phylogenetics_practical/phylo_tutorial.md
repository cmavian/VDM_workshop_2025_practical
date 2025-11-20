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


8. Examples for running tree

```

iqtree -nt 2 -s input.mafft.fst -m TEST -mrate R -B 1000

```

Explanation of tags:

- `-nt` determines the number of threads

- `-s` defines your alignment file *(required)*

- `-m` defines the substitution model
	
	- if you use `TEST` it will test what the best model is automatically
	
- `-mrate R` *(this is not necessary to know)* historically a gamma distribution
has been used for assessing rate heterogeneity across the alignment. However, recent
research has shown that gamma distributions introduce small (or large) biases ([read here for details](https://www.biorxiv.org/content/10.1101/2024.08.01.606208v1.full-text))
Instead, free rates is a better option. Use `-mrate R` to restrict your
model search to free rate models.

- `-B 1000` number of ultrafast bootstrap iterations. If you want to use real
bootstraps try `-b` (this is much much much slower)

