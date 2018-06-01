from ppl import *
import operator, math, heapq

class mVSA:
    
    def __init__(self, vectors):
        assert len(vectors) > 0, "There must be at least 1 vector"
        vec_len = [len(vec) for vec in vectors]
        assert vec_len[1:] == vec_len[:-1], "All vectors must have same length"
        self.dim = vec_len[0]
        for vec in vectors:
            assert vec != [0] * self.dim, "All vectors must be nonzero" 
        self.vectors = vectors
        self.axis = [Variable(i) for i in range(0, self.dim)]
        self.hyperplanes = None
        self.eps = (1e-8, 1e-3)
        vec_num = len(vectors)
        self.sub_vectors = [[vectors[i][k]-vectors[j][k] for k in range(0, self.dim)] 
                            for i in range(0, vec_num) for j in range(i+1, vec_num)]
        self.sub_vectors = filter(lambda x : x != [0] * self.dim, self.sub_vectors)
        self.SMA_reps = self.SMA_representatives(self.sub_vectors)
    
    def norm(self, vec):
        return math.sqrt(sum([vec[i]*vec[i] for i in range(0, self.dim)]))
        
    def sum_norm(self, vectors):
        vec_sum = [sum([vec[i] for vec in vectors]) for i in range(0, self.dim)]
        return self.norm(vec_sum)
        
    def hyperplane(self, _vec):
        vec = [v for v in _vec]
        divisor = [1 for v in vec]
        for i in range(0, self.dim):
            while abs(round(vec[i]) - vec[i]) > > self.eps[0] * divisor[i]:
                vec[i] *= 10
                divisor[i] *= 10
            vec[i] = round(vec[i])
        max_divisor = max(divisor)
        return sum([(vec[i]*max_divisor/divisor[i]) * self.axis[i] for i in range(0, self.dim)])
    
    def linear_order(self, vec):
        def linear_order_v(vec1, vec2):
            v = sum([(vec1[i]-vec2[i])*vec[i] for i in range(0, self.dim)])
            return 1 if v > 0 else -1 if v < 0 else 0
        return linear_order_v
    
    def support_point(self, generators):
        for generator in generators:
            if generator.is_point() and generator.space_dimension() != 0:
                supp_p = [(coef+0.0)/generator.divisor() for coef in generator.coefficients()]
                vec_norm = self.norm(supp_p)
                return tuple(float(coef/vec_norm) for coef in supp_p)
        return None
    
    def intersect(self, generators, vec):
        for generator in generators:
            if generator.is_ray():
                ray_coord = generator.coefficients()
                inner_product = sum([ray_coord[i]*vec[i] for i in range(0, self.dim)]) + 0.0
                if abs(inner_product)/(self.norm(vec)*self.norm(ray_coord)) < self.eps[1]:
                    return True
        return False
    
    def SMA_representative(self, idx_vec, v_idx, cs, hps, SMA_reps):
        if idx_vec[0] == v_idx:
        	v_idx += 1
            hps.append(True)
        if v_idx >= len(self.hyperplanes):
            generators = NNC_Polyhedron(cs).minimized_generators()
            SMA_reps.append((self.support_point(generators), hps))
            cs_inv = Constraint_System()
            for c in cs:
                c_inv = sum([c.coefficient(self.axis[i])*self.axis[i] for i in range(0, self.dim)]) < 0
                cs_inv.insert(c_inv)
            generators_inv = NNC_Polyhedron(cs_inv).minimized_generators()
            hps_inv = [False if hp else True for hp in hps]
            SMA_reps.append((self.support_point(generators_inv), hps_inv))
        else:
            cs_pos = Constraint_System(cs)
            cs_neg = Constraint_System(cs)
            cs_pos.insert(self.hyperplanes[v_idx] > 0)
            cs_neg.insert(self.hyperplanes[v_idx] < 0)
            generators_pos = NNC_Polyhedron(cs_pos).minimized_generators()
            generators_neg = NNC_Polyhedron(cs_neg).minimized_generators()
            if self.intersect(generators_pos, idx_vec[1]):
                self.SMA_representative(idx_vec, v_idx+1, cs_pos, hps+[True], SMA_reps)
            if self.intersect(generators_neg, idx_vec[1]):
                self.SMA_representative(idx_vec, v_idx+1, cs_neg, hps+[False], SMA_reps)

    def SMA_representatives(self, vectors):
        self.hyperplanes, SMA_reps = [], []
        for vec in vectors:
            self.hyperplanes.append(self.hyperplane(vec))
        for idx, vec in enumerate(vectors):
            self.SMA_representative((idx, vec), 0, self.hyperplanes[idx] > 0, [], SMA_reps) 
        SMA_reps_uniq = []
        for i in range(0, len(SMA_reps)):
            duplicate = False
            for j in range(i+1, len(SMA_reps)):
                if SMA_reps[i][1] == SMA_reps[j][1]:
                    duplicate = True
                    break
            if not duplicate:
                SMA_reps_uniq.append(SMA_reps[i][0])
        return SMA_reps_uniq
    
    def solve(self, M, average, top_k = 1, inverse = False):
        vec_num = len(self.vectors)
        if vec_num == 1:
            return [(self.sum_norm(self.vectors), self.vectors)]
        if top_k > 1:
            assert len(self.sub_vectors)*2 == vec_num*(vec_num-1), "All vectors must be unique when top_k > 1"
        candidates = []
        for SMA_rep in self.SMA_reps:
            indexed_vectors = zip(self.vectors, range(0, vec_num))
            indexed_vectors.sort(self.linear_order(SMA_rep), key = operator.itemgetter(0), reverse = not inverse)
            for m in M:
                [vectors, indices] = zip(*indexed_vectors[0:m])
                vec_sum_norm = self.sum_norm(vectors)
                if average:
                    vec_sum_norm /= float(len(vectors))
                candidates.append((vec_sum_norm, sorted(indices)))
        candidates.sort(key = operator.itemgetter(1))
        for i in range(0, len(candidates)-1):
            if candidates[i][1] == candidates[i+1][1]:
                candidates[i] = None
        candidates = filter(lambda x: x is not None, candidates)
        top_k = min(top_k, len(candidates))
        if inverse:
            return heapq.nsmallest(top_k, candidates, key = operator.itemgetter(0))
        else:
            return heapq.nlargest(top_k, candidates, key = operator.itemgetter(0))
        
    def m_VS(self, m, top_k = 1, inverse = False):
        assert m <= len(self.vectors), "m must not be greater than #vectors"
        return self.solve([m], False, top_k, inverse)
        
    def VSA(self, top_k = 1, inverse = False):
        return self.solve(range(1, len(self.vectors)+1), True, top_k, inverse)

