//
//  mVSA.hpp
//  Solution implementation for m-VS and VSA problems
//
//  Created by Ji, Yi on 21/05/2018.
//  Copyright © 2018 jiyi. All rights reserved.
//

#include "ppl.hh"
#include <vector>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream> // remove this after debug

using namespace Parma_Polyhedra_Library;

template <typename T>
class mVSA
{
    
public:
    
    template <typename V>
    struct value_trait
    {
        typedef typename V::value_type value_type;
    };
    
    template <typename V>
    struct value_trait<V*>
    {
        typedef V value_type;
    };
    
    template <typename V>
    struct value_trait<V**>
    {
        typedef V* value_type;
    };
    
    typedef typename std::remove_const<typename value_trait<typename value_trait<T>::value_type>::value_type>::type value_t;
    typedef typename std::vector<std::shared_ptr<const value_t>>::const_iterator iter_t;
    typedef typename std::vector<std::shared_ptr<const value_t>> vec_ptr_t;
    typedef unsigned int idx_t;
    
    mVSA(const T& arg_vectors, const idx_t arg_dim) : dim(arg_dim)
    {
        for (typename T::const_iterator vi = arg_vectors.begin(); vi != arg_vectors.end(); ++vi)
        {
            std::shared_ptr<value_t> vec(new value_t[dim], std::default_delete<value_t[]>());
            for (idx_t i = 0; i < dim; ++i)
                vec.get()[i] = (*vi)[i];
            vectors.push_back(vec);
        }
        init();
    }
    
    mVSA(const T& arg_vectors, const idx_t vec_num, const idx_t arg_dim) : dim(arg_dim)
    {
        for (idx_t i = 0; i < vec_num; ++i)
        {
            std::shared_ptr<value_t> vec(new value_t[dim], std::default_delete<value_t[]>());
            for (idx_t j = 0; j < dim; ++j)
                vec.get()[j] = arg_vectors[i][j];
            vectors.push_back(vec);
        }
        init();
    }
    
    inline const mpf_class inner_product(std::shared_ptr<const mpf_class> vec1,
                                         std::shared_ptr<const mpf_class> vec2) const
    {
        mpf_class inner_product = 0;
        for (idx_t i = 0; i < dim; ++i)
            inner_product += vec1.get()[i] * vec2.get()[i];
        return inner_product;
    }
    
    inline const mpf_class norm(std::shared_ptr<const mpf_class> vec) const
    {
        return sqrt(inner_product(vec, vec));
    }
    
    const mpf_class sum_norm(const vec_ptr_t& arg_vectors) const
    {
        std::shared_ptr<mpf_class> vec(new mpf_class[dim], std::default_delete<value_t[]>());
        for (idx_t i = 0; i < dim; ++i)
            vec.get()[i] = mpf_class(0);
        for (iter_t vi = arg_vectors.begin(); vi != arg_vectors.end(); ++vi)
        {
            for (idx_t i = 0; i < dim; ++i)
                vec.get()[i] += mpf_class(vi->get()[i]);
        }
        return norm(vec);
    }
    
    Linear_Expression hyperplane(const std::shared_ptr<const value_t> _vec) const
    {
        std::shared_ptr<value_t> vec(new value_t[dim], std::default_delete<value_t[]>());
        std::shared_ptr<mpz_class> divisor(new mpz_class[dim], std::default_delete<mpz_class[]>());
        for (idx_t i = 0; i < dim; ++i)
        {
            vec.get()[i] = _vec.get()[i];
            divisor.get()[i] = 1;
        }
        mpz_class max_divisor = 1;
        for (idx_t i = 0; i < dim; ++i)
        {
            while (abs(round(vec.get()[i]) - vec.get()[i]) > eps[0] * divisor.get()[i])
            {
                vec.get()[i] *= 10;
                divisor.get()[i] *= 10;
            }
            vec.get()[i] = round(vec.get()[i]);
            max_divisor = std::max(divisor.get()[i], max_divisor);
        }
        Linear_Expression le;
        for (idx_t i = 0; i < dim; ++i)
            le += (vec.get()[i] * max_divisor / divisor.get()[i]) * axis[i];
        le.print();
        std::cout << std::endl;
        return le;
    }
    
    std::shared_ptr<mpz_class> coord(const Generator& g) const
    {
        std::shared_ptr<mpz_class> g_coord(new mpz_class[dim], std::default_delete<mpz_class[]>());
        for (idx_t i = 0; i < dim; ++i)
            g_coord.get()[i] = mpz_class(g.coefficient(axis[i]));
        return g_coord;
    }
    
    std::shared_ptr<mpz_class> support_point(const Generator_System& gs) const
    {
        for (Generator_System::const_iterator gi = gs.begin(); gi != gs.end(); ++gi)
        {
            if (gi->is_point() && gi->space_dimension() != 0)
                return coord(*gi);
        }
        return nullptr;
    }
    
    bool intersect(const Generator_System& gs, const std::shared_ptr<const value_t> vec) const
    {
        for (Generator_System::const_iterator gi = gs.begin(); gi != gs.end(); ++gi)
        {
            if (gi->is_ray())
            {
                std::shared_ptr<mpf_class> ray(new mpf_class[dim], std::default_delete<mpf_class[]>());
                std::shared_ptr<mpf_class> vec_f(new mpf_class[dim], std::default_delete<mpf_class[]>());
                for (idx_t i = 0; i < dim; ++i)
                {
                    ray.get()[i] = mpf_class(gi->coefficient(axis[i]));
                    vec_f.get()[i] = mpf_class(vec.get()[i]);
                }
                const mpf_class ray_vec_product = inner_product(ray, vec_f);
                if (abs(ray_vec_product)/(norm(ray)*norm(vec_f)) < eps[1])
                    return true;
            }
        }
        return false;
    }
    
    void SMA_representative(const vec_ptr_t& arg_vectors, idx_t vec_idx, idx_t v_idx, Constraint_System cs, std::vector<bool> hps,
                            std::vector<std::pair<std::shared_ptr<mpz_class>, std::vector<bool>>>& arg_SMA_reps) const
    {
        if (v_idx == vec_idx)
        {
            v_idx++;
            hps.push_back(true);
        }
        if (v_idx >= hyperplanes.size())
        {
            Generator_System gs = NNC_Polyhedron(cs).minimized_generators();
            Constraint_System cs_inv;
            for (Constraint_System::const_iterator ci = cs.begin(); ci != cs.end(); ++ci)
            {
                Linear_Expression c_inv;
                for (idx_t i = 0; i < dim; ++i)
                    c_inv += ci->coefficient(axis[i]) * axis[i];
                cs_inv.insert(c_inv < 0);
            }
            Generator_System gs_inv = NNC_Polyhedron(cs_inv).minimized_generators();
            std::vector<bool> hps_inv;
            for (std::vector<bool>::iterator i = hps.begin(); i != hps.end(); ++i)
                hps_inv.push_back(!*i);
            std::pair<std::shared_ptr<mpz_class>, std::vector<bool>> rep(support_point(gs), hps),
                                                                     rep_inv(support_point(gs_inv), hps_inv);
            arg_SMA_reps.push_back(rep);
            arg_SMA_reps.push_back(rep_inv);
        }
        else
        {
            Constraint_System cs_pos(cs), cs_neg(cs);
            cs_pos.insert(hyperplanes[v_idx] > 0);
            cs_neg.insert(hyperplanes[v_idx] < 0);
            Generator_System gs_pos = NNC_Polyhedron(cs_pos).minimized_generators();
            Generator_System gs_neg = NNC_Polyhedron(cs_neg).minimized_generators();
            if (intersect(gs_pos, arg_vectors[vec_idx]))
            {
                std::vector<bool> hps_pos(hps);
                hps_pos.push_back(true);
                SMA_representative(arg_vectors, vec_idx, v_idx+1, cs_pos, hps_pos, arg_SMA_reps);
            }
            if (intersect(gs_neg, arg_vectors[vec_idx]))
            {
                std::vector<bool> hps_neg(hps);
                hps_neg.push_back(false);
                SMA_representative(arg_vectors, vec_idx, v_idx+1, cs_neg, hps_neg, arg_SMA_reps);
            }
        }
    }
    
    std::vector<std::shared_ptr<mpz_class>> SMA_representatives(const vec_ptr_t& arg_vectors)
    {
        std::vector<std::pair<std::shared_ptr<mpz_class>, std::vector<bool>>> arg_SMA_reps;
        hyperplanes.clear();
        for (iter_t vi = arg_vectors.begin(); vi != arg_vectors.end(); ++vi)
            hyperplanes.push_back(hyperplane(*vi));
        for (idx_t i = 0; i < dim; ++i)
            SMA_representative(arg_vectors, i, 0, Constraint_System(hyperplanes[i] > 0), std::vector<bool>(), arg_SMA_reps);
        //print_SMA_reps(arg_SMA_reps);
        std::vector<std::shared_ptr<mpz_class>> SMA_reps_uniq;
        for (auto i = arg_SMA_reps.begin(); i != arg_SMA_reps.end(); ++i)
        {
            auto j = i; ++j;
            for (; j != arg_SMA_reps.end(); ++j)
            {
                bool same_cone = true;
                for (idx_t k = 0; k < hyperplanes.size() && same_cone; ++k)
                    same_cone = same_cone && (i->second[k] == j->second[k]);
                if (same_cone)
                    break;
            }
            if (j == arg_SMA_reps.end())
                SMA_reps_uniq.push_back(i->first);
        }
        print_SMA_reps(SMA_reps_uniq);
        return SMA_reps_uniq;
    }
    
    void print_SMA_reps(const std::vector<std::shared_ptr<mpz_class>>& arg_SMA_reps)
    {
        std::cout << arg_SMA_reps.size() << std::endl;
        for (idx_t i = 0; i < arg_SMA_reps.size(); ++i)
        {
            std::cout << "[";
            for (idx_t j = 0; j < dim; ++j)
                std::cout << arg_SMA_reps[i].get()[j].get_str() << ", ";
            std::cout << "]" << std::endl;
        }
    }
    
    void print_SMA_reps(const std::vector<std::pair<std::shared_ptr<mpz_class>, std::vector<bool>>>& arg_SMA_reps)
    {
        std::cout << arg_SMA_reps.size() << std::endl;
        for (idx_t i = 0; i < arg_SMA_reps.size(); ++i)
        {
            std::cout << "[";
            for (idx_t j = 0; j < dim; ++j)
                std::cout << arg_SMA_reps[i].first.get()[j].get_str() << ", ";
            std::cout << "], ";
            std::vector<bool> hps = arg_SMA_reps[i].second;
            for (std::vector<bool>::iterator vi = hps.begin(); vi != hps.end(); ++vi)
            {
                std::cout << *vi << ", ";
            }
            std::cout << std::endl;
        }
    }
    
    void init()
    {
        assert(vectors.size() > 0 && dim > 0);
        for (iter_t vi = vectors.begin(); vi != vectors.end(); ++vi)
        {
            bool nonzero = false;
            for (idx_t i = 0; i < dim; ++i)
                nonzero = nonzero || vi->get()[i] != 0;
            assert(nonzero);
            for (iter_t vj = vi+1; vj != vectors.end(); ++vj)
            {
                std::shared_ptr<value_t> sub_vec(new value_t[dim], std::default_delete<value_t[]>());
                bool unique = false;
                for (idx_t i = 0; i < dim; ++i)
                {
                    sub_vec.get()[i] = vi->get()[i] - vj->get()[i];
                    unique = unique || sub_vec.get()[i] != 0;
                }
                if (unique)
                    sub_vectors.push_back(sub_vec);
            }
        }
        for (idx_t i = 0; i < dim; ++i)
            axis.push_back(Variable(i));
        SMA_reps = SMA_representatives(sub_vectors);
    }
    
private:
    
    vec_ptr_t vectors, sub_vectors;
    std::vector<std::shared_ptr<mpz_class>> SMA_reps;
    std::vector<Variable> axis;
    std::vector<Linear_Expression> hyperplanes;
    const idx_t dim;
    const mpf_class eps[2] = {1e-6, 1e-3};
};
