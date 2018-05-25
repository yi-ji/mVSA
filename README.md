# mVSA

This repo contains implementation of solutions for some vector subset choice (or vector sum) problems, namely:

- `m-VS`: Vector Subset with the Maximum/Minimum Sum Norm (*m* stands for given number of vectors, optional) <br/>
- `VSA`: Vector Subset with the Maximum/Minumum Averaged Square of the Sum Norm

also enriching the original algorithms with top-*k* extension.

## Python Implementation

### Requirement
- [pplpy](https://github.com/videlec/pplpy)

### Usage

```python
>>> import mVSA

>>> vectors = [[1, 0, 5.1], [5.2, 0, 0], [1.1, 2, 3], [2, 5, 4]] # supporting both integer & float
>>> mvsa = mVSA.mVSA(vectors)

>>> print mvsa.m_VS(m = 2, top_k = 3)
[(10.807867504739313, [0, 3]), (10.373523991392704, [2, 3]), (9.635351576356724, [1, 3])]
>>> print mvsa.m_VS(m = 3, top_k = 2, inverse = True) # inverse for minimum case
[(11.086027241532468, [0, 1, 2]), (12.918591254467339, [1, 2, 3])]
>>> print mvsa.VSA(top_k = 2, inverse = True)
[(3.6293938887919017, [1, 2]), (3.695342413844156, [0, 1, 2])]
>>> print mvsa.VSA(top_k = 4)
[(6.708203932499369, [3]), (5.403933752369657, [0, 3]), (5.2, [1]), (5.197114584074513, [0])]

```

## C++ Implementation

TBC

## References

- [On Polynomial Solvability of Some Problems of a Vector Subset Choice in a Euclidean Space of Fixed Dimension](https://link.springer.com/article/10.1134/S1990478910010084)
- [Polynomial Algorithms for Solving the Vector Sum Problem](https://link.springer.com/article/10.1134/S1990478907030027)
