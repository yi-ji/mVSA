//
//  test_mVSA.cpp
//  Vector-Sum-Problems
//
//  Created by Ji, Yi on 02/06/2018.
//  Copyright © 2018 jiyi. All rights reserved.
//

#include <array>
#include <random>
#include "mVSA.hpp"

using namespace std;

void test_equal(const vector<pair<double, vector<unsigned int>>>& ans1,
                const vector<pair<double, vector<unsigned int>>>& ans2)
{
    assert(ans1.size() == ans2.size());
    vector<pair<double, vector<unsigned int>>>::const_iterator i, j;
    for (i = ans1.begin(), j = ans2.begin(); i != ans1.end() && j != ans2.end(); ++i, ++j)
    {
        assert(abs(i->first-j->first) < 1e-3);
        assert(i->second.size() == j->second.size());
    }
}

int main()
{
    const float vs[5][4] = {{0.1, 2.3, 4.2, 4.6}, {12, 3.4, 4.33, 5.555}, {1, 3, 5, 6}, {4.2, 5, 9, 8}, {1, 4, 5.5, 6.4}};
 
    array<const float*, 5> vectors = {vs[0], vs[1], vs[2], vs[3], vs[4]};
    mVSA<array<const float*, 5>> mvsa(vectors, 4);

    test_equal(mvsa.m_VS(3, 3, true), mvsa.m_VS_brute_force(3, 3, true));
    test_equal(mvsa.VSA(5), mvsa.VSA_brute_force(5));
    
    // using ramdomly generatored vectors and parameters to test
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dim(2, 8), num(2, 16), top_k(1, 4), value(-100, 100);
    for (int t = 0; t < 100; ++t)
    {
        vector<vector<int>> random_vectors;
        int vec_num = num(rng);
        int vec_dim = dim(rng);
        for (int i = 0; i < vec_num; ++i)
        {
            vector<int> vec;
            for (int j = 0; j < vec_dim; ++j)
                vec.push_back(value(rng));
            random_vectors.push_back(vec);
        }
        std::uniform_int_distribution<int> M(1, vec_num);
        mVSA<vector<vector<int>>> r_mvsa(random_vectors, vec_dim);
        int m = M(rng), k = top_k(rng);
        bool inverse = value(rng) % 2 == 0;
        test_equal(r_mvsa.m_VS(m, k, inverse), r_mvsa.m_VS_brute_force(m, k, inverse));
        test_equal(r_mvsa.VSA(k, inverse), r_mvsa.VSA_brute_force(k, inverse));
    }
    
    cout << "All tests passed" << endl;
}
