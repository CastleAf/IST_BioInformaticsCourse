# Smith-Waterman Algorithm

**Afonso Castelão, ist190700, MEIC-A**

Bioinformatics course's first lab delivery.

I was assign to implement the Smith-Waterman algorithm so I could align two pairs of amino acid sequences. Algorithm is fully working with BLOSUM50 as the score matrix.

#### How to run:
```
$ python3 smith-waterman.py
```
Simply run the command above and input two strings sequences and a negative integer.

For example:
```
$ python3 smith-waterman.py
First Amino Acid Sequence: WEXIWEW
Second Amino Acid Sequence: PWEWWEW
Gap Penalty Value: -8
```

---

#### Here's a graphical representation of the example above:
![Smith-Waterman alignment Example](https://github.com/CastleAf/IST_BioInformaticsCourse/blob/master/Smith-Waterman%20Algorithm/img/Example_Alignment.png?raw=true | width=48)

Below there's a link that does graphical representation of the algorithm. When using it, make sure you select Smith-Waterman alignment type (this link's program is not mine):

* https://gtuckerkellogg.github.io/pairwise/demo/

