# mVSA

This repo contains implementation of solutions for some vector subset choice (or vector sum) problems, namely:

- `m-VS`: *Vector Subset with the Maximum/Minimum Sum Norm* (`m` stands for given number of vectors, optional) <br/>
- `VSA`: *Vector Subset with the Maximum/Minimum Averaged Square of the Sum Norm*

also enriching the original algorithms with top-*k* extension.

All input vectors must be nonzero and unique (when *k* > 1).

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

### Requirement
- [Parma Polyhedra Library](http://www.bugseng.com/ppl-download) (PPL)

### Usage

Instantiation by any template parameter within `mVSA<{VectorSet<Subscriptable<T>>`, `VectorSet<T*>`, `T**}>` can be accepted. </br>
`VectorSet` allows all C++ STL `SequenceContainer` and some `AssociativeContainer`. </br>
`Subscriptable` must support subscript [] operator for vector indexing. </br>
The vectors dimension must be provided as a constructor parameter.

Output result is an ordered `std::vector<std::pair<double, std::vector<unsigned int>>>` where each pair consists of the max/min sum norm value and corresponding vector indices.

```C++
#include <set>
#include <array>
#include <vector>
#include "mVSA.hpp"

using namespace std;

int main() // some examples
{
    const float vs[5][4] = {{0.1, 2.3, 4.2, 4.6}, {12, 3.4, 4.33, 5.555}, {1, 3, 5, 6}, {4.2, 5, 9, 8}, {1, 4, 5.5, 6.4}};
    const int vs2[4][3] = {{1, 0, 5}, {5, 0, 0}, {1, 2, 3}, {2, 5, 4}};
    double vs3[4][3] = {{1, 0, 5.1}, {5.2, 0, 0}, {1.1, 2, 3}, {2, 5, 4}};
    
    
    array<const float*, 5> vectors1 = {vs[0], vs[1], vs[2], vs[3], vs[4]};
    mVSA<array<const float*, 5>> mvsa1(vectors1, 4);
    mvsa1.print_result(mvsa1.m_VS(3, 3, true));
    // output: [(24.4129, [0, 2, 4, ]), (26.2934, [0, 1, 2, ]), (27.1405, [0, 1, 4, ]), ]
    mvsa1.print_result(mvsa1.VSA(5));
    // output: [(14.3236, [1, ]), (13.6982, [3, ]), (13.176, [1, 3, ]), (11.5584, [1, 3, 4, ]), (11.4635, [3, 4, ]), ]
    
    
    vector<vector<int>> vectors2 = {{1, 0, 5}, {5, 0, 0}, {1, 2, 3}, {2, 5, 4}};
    mVSA<vector<vector<int>>> mvsa2(vectors2, 3);
    mvsa2.print_result(mvsa2.m_VS(2, 3));
    // output: [(10.7238, [0, 3, ]), (10.3441, [2, 3, ]), (9.48683, [1, 3, ]), ]
    mvsa2.print_result(mvsa2.VSA(2, true));
    // output: [(3.5, [1, 2, ]), (3.60555, [0, 1, 2, ]), ]
    
    
    set<double*> vectors3 = {vs3[0], vs3[1], vs3[2], vs3[3]};
    mVSA<set<double*>> mvsa3(vectors3, 3);
    mvsa3.print_result(mvsa3.VSA(4));
    // output: [(6.7082, [3, ]), (5.40393, [0, 3, ]), (5.2, [1, ]), (5.19711, [0, ]), ]
    
    const int* vsp[] = {vs2[0], vs2[1], vs2[2], vs2[3]};
    mVSA<const int**> mvsa4(vsp, 4, 3);
    mvsa4.print_result(mvsa4.m_VS(3, 2, true));
    // output: [(10.8167, [0, 1, 2, ]), (12.7279, [1, 2, 3, ]), ]
    mvsa4.print_result(mvsa4.VSA(4));
    // output: [(6.7082, [3, ]), (5.3619, [0, 3, ]), (5.17204, [2, 3, ]), (5.09902, [0, ]), ]
    
    return 0;
}
```

## References

- [On Polynomial Solvability of Some Problems of a Vector Subset Choice in a Euclidean Space of Fixed Dimension](https://link.springer.com/article/10.1134/S1990478910010084)
- [Polynomial Algorithms for Solving the Vector Sum Problem](https://link.springer.com/article/10.1134/S1990478907030027)
