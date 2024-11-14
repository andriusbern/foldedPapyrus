**Folded Papyrus** 
This dataset contains [AlphaFold](https://www.nature.com/articles/s41586-021-03819-2)'s predictions and internal latent representations for a large subset* of the proteins that are referenced in [Papyrus][1] (a large scale curated dataset aimed at bioactivity predictions). 


### **Dataset structure**
---
<!-- <sub> -->
<!-- <sup> -->

&#8595; **data/*** -> main data directory

> &#8595; **data/PID/***         -> data of a single protein of length **L**
>
> > 
> > | Filename | Description |Tensor shape| Lightweight |
> > | --- | --- | --- | --- |
> > | **single.npy***      | ($s_i$)  evoformer single representation   | *[**L** x 384]* |   ✔️
> > | **structure.npy***   | $(a_i)$  output of the last layer of structure module  | *[**L** x 384]* | ✔️
> > | **msa.npy***         | $(m_{si})$ processed MSA representation            |  *[**N** x **L** x 256]* |
> > | **pair.npy***        |  $(z_{ij}$) evoformer pair representation  |  *[**L** x **L** x 128]* | 
> > | **PID.pdb**     | 3D protein structure prediction |  | ✔️
> > | **PID_unrelaxed.pdb**     | 3D protein structure prediction w/o relaxation step (D) |  | ✔️
> > | **confidence.npy***   | confidence in structure prediction (0-100)  | *[1]* | ✔️
> > | **plldt.npy***   | confidence in structure prediction per residue  | *[**L**]* | ✔️
> > |**PID.fasta**    | protein amino acid sequence and metadata | |✔️ 
> > | **timings.json**    | Processing log  | | ✔️
<!-- > > | **distogram_logits.npy***        | probability distribution of the structure distogram  |  *[**L** x **L** x 64]* |  -->
<!-- > > | **features.pkl**    | Additional metrics that are saved automatically | | ✔️    -->
> &#8595; **data/PID2/***  -> data of protein #2
> > **...**
> > 
<!-- </sup> -->
<!-- </sub> -->

**Note: **L**: sequence length, **N**: number of aligned sequences via MSA.*

---

<br /> 
---

## References:
1. Béquignon OJM, Bongers BJ, Jespers W, IJzerman AP, van de Water B, van Westen GJP. [Papyrus - A large scale curated dataset aimed at bioactivity predictions. ](https://chemrxiv.org/engage/chemrxiv/article-details/617aa2467a002162403d71f0)ChemRxiv. Cambridge: Cambridge Open Engage; 2021; This content is a preprint and has not been peer-reviewed.
2. Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

[1]: <https://chemrxiv.org/engage/chemrxiv/article-details/617aa2467a002162403d71f0> "Béquignon OJM, Bongers BJ, Jespers W, IJzerman AP, van de Water B, van Westen GJP. Papyrus - A large scale curated dataset aimed at bioactivity predictions. ChemRxiv. Cambridge: Cambridge Open Engage; 2021; This content is a preprint and has not been peer-reviewed."
