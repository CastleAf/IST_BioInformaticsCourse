# HMM Algorithms

**Afonso Castelão, ist190700, MEIC-A**

Bioinformatics course's second lab delivery.

I was assign to implement some Hidden Markov Models Algorithms so I could solve some presented problems. For this report, I coded three algorithms: **Viterbi Algorithm, Forward Algorithm and Backward Algorithm**. Finally, I used the Posterior Probability formula to solve the posterior decoding exercise.

#### Here's the model I've worked with in this program:
![HMM Transition Model](https://github.com/CastleAf/IST_BioInformaticsCourse/blob/master/HMM_Algorithms/img/HMM_Transitions.svg?raw=true)


**To Note:** As you can observe, I used the same Hidden Markov Model for all the functions. This program isn't really that useful because it only works for this given HMM. Therefore, this code isn't reusable, it only works as an example of how the algorithms work.

---

#### How to run:
```
$ python3 HMM_Algorithms.py
```
Simply run the command above and input a sequence of DNA symbols observations.

For example:
```
$ python3 HMM_Algorithms.py
Please input the DNA sequence: ACGTA
```


