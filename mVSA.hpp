//
//  mVSA.hpp
//  Solution implementation for m-VS and VSA problems
//
//  Created by Ji, Yi on 21/05/2018.
//  Copyright © 2018 jiyi. All rights reserved.
//

#ifndef mVSA_HPP
#define mVSA_HPP

#include "ppl.hh"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream> // for debug

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
    typedef typename std::vector<std::pair<double, std::vector<idx_t>>> result_t;
    
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
    
    const double sum_norm(const vec_ptr_t& arg_vectors) const
    {
        std::shared_ptr<mpf_class> vec(new mpf_class[dim], std::default_delete<mpf_class[]>());
        for (idx_t i = 0; i < dim; ++i)
            vec.get()[i] = mpf_class(0);
        for (iter_t vi = arg_vectors.begin(); vi != arg_vectors.end(); ++vi)
        {
            for (idx_t i = 0; i < dim; ++i)
                vec.get()[i] += mpf_class(vi->get()[i]);
        }
        return norm(vec).get_d();
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
                if (abs(ray_vec_product) / (norm(ray)*norm(vec_f)) < eps[1])
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
        for (idx_t i = 0; i < arg_vectors.size(); ++i)
            SMA_representative(arg_vectors, i, 0, Constraint_System(hyperplanes[i] > 0), std::vector<bool>(), arg_SMA_reps);
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
        return SMA_reps_uniq;
    }
    
    // for debug
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
    
    // for debug
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
                std::cout << *vi << ", ";
            std::cout << std::endl;
        }
    }
    
    // for debug
    void print_result(const result_t& result)
    {
        std::cout << "[";
        for (typename result_t::const_iterator i = result.begin(); i != result.end(); ++i)
        {
            std::cout << "(" << i->first << ", [";
            for (std::vector<idx_t>::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
                std::cout << *j << ", ";
            std::cout << "]), ";
        }
        std::cout << "]" << std::endl;
    }
    
    void init()
    {
        assert(vectors.size() > 0 && dim > 0 && "There must be at lease 1 vector");
        for (iter_t vi = vectors.begin(); vi != vectors.end(); ++vi)
        {
            bool nonzero = false;
            for (idx_t i = 0; i < dim; ++i)
                nonzero = nonzero || vi->get()[i] != 0;
            assert(nonzero && "All vectors must be nonzero");
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
    
    std::function<bool(std::pair<double, std::vector<idx_t>>, std::pair<double, std::vector<idx_t>>)> norm_order(bool inverse)
    {
        return [inverse] (auto v1, auto v2) -> bool
        {
            bool cmp = v1.first < v2.first;
            return inverse ? !cmp : cmp;
        };
    }
    
    std::function<bool(std::pair<double, std::vector<idx_t>>, std::pair<double, std::vector<idx_t>>)> space_order()
    {
        return [] (auto v1, auto v2) -> bool
        {
            if (v1.second.size() != v2.second.size())
                return v1.second.size() < v2.second.size() ? true : false;
            for (idx_t i = 0; i < v1.second.size(); ++i)
            {
                if (v1.second[i] < v2.second[i])
                    return true;
                else if (v1.second[i] > v2.second[i])
                    return false;
            }
            return false;
        };
    }
    
    std::function<bool(std::pair<std::shared_ptr<const value_t>, idx_t>, std::pair<std::shared_ptr<const value_t>, idx_t>)> linear_order(const std::shared_ptr<mpz_class> SMA_rep)
    {
        return [this, SMA_rep] (auto v1, auto v2) -> bool
        {
            mpf_class order_v(0);
            for (idx_t i = 0; i < dim; ++i)
                order_v += mpf_class(v1.first.get()[i] - v2.first.get()[i]) * mpf_class(SMA_rep.get()[i]);
            return order_v > 0 ? true : false;
        };
    }
    
    result_t select_top_k(result_t& result, idx_t top_k, bool inverse)
    {
        result_t top_k_result;
        top_k = std::min(top_k, (idx_t) result.size());
        auto heap_cmp = norm_order(inverse);
        std::make_heap(result.begin(), result.end(), heap_cmp);
        for (idx_t k = 0; k < top_k; ++k)
        {
            top_k_result.push_back(result.front());
            std::pop_heap(result.begin(), result.end(), heap_cmp);
            result.pop_back();
        }
        return top_k_result;
    }
    
    result_t solve(std::vector<idx_t> M, bool average, idx_t top_k, bool inverse)
    {
        if (vectors.size() == 1)
            return result_t(1, std::pair<double, std::vector<idx_t>>(sum_norm(vectors), std::vector<idx_t>(1, 0)));
        if (top_k > 1)
            assert(sub_vectors.size()*2 == vectors.size()*(vectors.size()-1) && "All vectors must be unique when top_k > 1");
        result_t candidates, candidates_uniq;
        for (std::vector<std::shared_ptr<mpz_class>>::iterator rep_i = SMA_reps.begin(); rep_i != SMA_reps.end(); ++rep_i)
        {
            std::vector<std::pair<std::shared_ptr<const value_t>, idx_t>> indexed_vectors;
            for (idx_t v = 0; v < vectors.size(); ++v)
                indexed_vectors.push_back(std::pair<std::shared_ptr<const value_t>, idx_t>(vectors[v], v));
            std::sort(indexed_vectors.begin(), indexed_vectors.end(), linear_order(*rep_i));
            for (std::vector<idx_t>::iterator m = M.begin(); m != M.end(); ++m)
            {
                std::vector<std::shared_ptr<const value_t>> m_vectors;
                std::vector<idx_t> indices;
                for (idx_t i = 0; i < *m; ++i)
                {
                    m_vectors.push_back(indexed_vectors[i].first);
                    indices.push_back(indexed_vectors[i].second);
                }
                std::sort(indices.begin(), indices.end());
                mpf_class vec_sum_norm(sum_norm(m_vectors));
                if (average)
                    vec_sum_norm /= *m;
                candidates.push_back(std::pair<double, std::vector<idx_t>>(vec_sum_norm.get_d(), indices));
            }
        }
        std::sort(candidates.begin(), candidates.end(), space_order());
        result_t::iterator i = candidates.begin(), j;
        candidates_uniq.push_back(*(i++));
        for (j = i-1; i != candidates.end(); j = i, ++i)
        {
            bool same_indices = i->second.size() == j->second.size();
            for (idx_t k = 0; k != i->second.size() && same_indices; ++k)
                same_indices = same_indices && (i->second[k] == j->second[k]);
            if (!same_indices)
                candidates_uniq.push_back(*i);
        }
        return select_top_k(candidates_uniq, top_k, inverse);
    }
    
    void m_power_set(idx_t m, idx_t idx, std::vector<idx_t> ps, std::vector<std::vector<idx_t>>& mps)
    {
        if (idx >= vectors.size() || ps.size() >= m)
        {
            if (ps.size() == m)
                mps.push_back(ps);
        }
        else
        {
            m_power_set(m, idx+1, ps, mps);
            ps.push_back(idx);
            m_power_set(m, idx+1, ps, mps);
        }
    }
    
    result_t m_VS_brute_force(idx_t m, idx_t top_k = 1, bool inverse = false)
    {
        result_t result;
        std::vector<std::vector<idx_t>> all_indices;
        m_power_set(m, 0, std::vector<idx_t>(), all_indices);
        for (std::vector<std::vector<idx_t>>::iterator i = all_indices.begin(); i != all_indices.end(); ++i)
        {
            vec_ptr_t arg_vectors;
            for (std::vector<idx_t>::iterator j = i->begin(); j != i->end(); ++j)
                arg_vectors.push_back(vectors[*j]);
            result.push_back(std::pair<double, std::vector<idx_t>>(sum_norm(arg_vectors), *i));
        }
        return select_top_k(result, top_k, inverse);
    }
    
    result_t VSA_brute_force(idx_t top_k = 1, bool inverse = false)
    {
        result_t result;
        for (idx_t m = 1; m <= vectors.size(); ++m)
        {
            result_t m_result = m_VS_brute_force(m, top_k, inverse);
            for (result_t::iterator i = m_result.begin(); i != m_result.end(); ++i)
            {
                i->first /= double(m);
                result.push_back(*i);
            }
        }
        return select_top_k(result, top_k, inverse);
    }
    
    result_t m_VS(idx_t m, idx_t top_k = 1, bool inverse = false)
    {
        assert(m <= vectors.size() && "m must not be greater than #vectors");
        return solve(std::vector<idx_t>(1, m), false, top_k, inverse);
    }
    
    result_t VSA(idx_t top_k = 1, bool inverse = false)
    {
        std::vector<idx_t> M;
        for (idx_t i = 1; i <= vectors.size(); ++i)
            M.push_back(i);
        return solve(M, true, top_k, inverse);
    }
    
private:
    
    vec_ptr_t vectors, sub_vectors;
    std::vector<std::shared_ptr<mpz_class>> SMA_reps;
    std::vector<Variable> axis;
    std::vector<Linear_Expression> hyperplanes;
    const idx_t dim;
    const mpf_class eps[2] = {1e-6, 1e-3};
};

#endif
