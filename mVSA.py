from ppl import *
import operator, math, heapq

class mVSA:
    
    def __init__(self, vectors):
        assert len(vectors) > 0, "There must be at least 1 vector"
        vec_len = [len(vec) for vec in vectors]
        assert vec_len[1:] == vec_len[:-1], "All vectors must have same length"
        self.vectors = vectors
        self.dim = vec_len[0]
        self.axis = [Variable(i) for i in range(0, self.dim)]
        vec_num = len(vectors)
        self.sub_vectors = [[vectors[i][k]-vectors[j][k] for k in range(0, self.dim)] 
                            for i in range(0, vec_num) for j in range(i+1, vec_num)]
        self.SMA_reps = self.SMA_representatives(self.sub_vectors)
        
    def is_close(self, a, b, epsilon = 1e-7):
        min_abs = min(abs(a), abs(b))
        return min_abs == 0 or abs(a-b)/(min_abs+0.0) < epsilon
        
    def sum_norm(self, vectors):
        vec_sum = [sum([vec[i] for vec in vectors]) for i in range(0, self.dim)]
        return math.sqrt(sum([vec_sum[i]*vec_sum[i] for i in range(0, self.dim)]))
    
    def to_integer(self, num):
        divisor = 1
        while not self.is_close(long(num), num):
            num *= 10
            divisor *= 10
        return divisor
        
    def hyperplane(self, vec):
        divisor = max([self.to_integer(vec[i]) for i in range(0, self.dim)])
        return sum([(vec[i]*divisor) * self.axis[i] for i in range(0, self.dim)])
    
    def linear_order(self, vec):
        def linear_order_v(vec1, vec2):
            v = sum([(vec1[i]-vec2[i])*vec[i] for i in range(0, self.dim)])
            return 1 if v > 0 else -1 if v < 0 else 0
        return linear_order_v
    
    def support_point(self, generators):
        for generator in generators:
            if generator.is_point() and generator.space_dimension() != 0:
                supp_p = [float((coef+0.0)/generator.divisor()) for coef in generator.coefficients()]
                norm = math.sqrt(sum([supp_p[i]*supp_p[i] for i in range(0, self.dim)]))
                return tuple(coef/norm for coef in supp_p)
        return None
    
    def intersect(self, generators, vec):
        for generator in generators:
            if generator.is_ray():
                ray_coord = generator.coefficients()
                if self.is_close(sum([ray_coord[i]*vec[i] for i in range(0, self.dim)]), 0):
                    return True
        return False
    
    def SMA_representative(self, vec, v_idx, vectors, cs):
        if v_idx < len(vectors) and vec == vectors[v_idx]:
            v_idx += 1
        if v_idx >= len(vectors):
            generators = NNC_Polyhedron(cs).minimized_generators()
            return [self.support_point(generators)]
        cs_pos = Constraint_System(cs)
        cs_neg = Constraint_System(cs)
        v_hp = self.hyperplane(vectors[v_idx])
        cs_pos.insert(v_hp > 0)
        cs_neg.insert(v_hp < 0)
        generators_pos = NNC_Polyhedron(cs_pos).minimized_generators()
        generators_neg = NNC_Polyhedron(cs_neg).minimized_generators()
        SMA_reps = []
        if self.intersect(generators_pos, vec):
            SMA_reps += self.SMA_representative(vec, v_idx+1, vectors, cs_pos)
        if self.intersect(generators_neg, vec):
            SMA_reps += self.SMA_representative(vec, v_idx+1, vectors, cs_neg)
        return SMA_reps

    def SMA_representatives(self, vectors):
        SMA_reps = []
        for vec in vectors:
            SMA_reps += self.SMA_representative(vec, 0, vectors, self.hyperplane(vec) > 0) 
        cs = Constraint_System()
        for vec in vectors:
            cs.insert(self.hyperplane(vec) < 0)
        generators = NNC_Polyhedron(cs).minimized_generators()
        supp_p = self.support_point(generators)
        if supp_p:
            SMA_reps.append(supp_p)
        SMA_reps_uniq = []
        for i in range(0, len(SMA_reps)):
            duplicate = False
            for j in range(i+1, len(SMA_reps)):
                if all([self.is_close(SMA_reps[i][k], SMA_reps[j][k]) for k in range(0, self.dim)]):
                    duplicate = True
                    break
            if not duplicate:
                SMA_reps_uniq.append(SMA_reps[i])
        return SMA_reps_uniq
    
    def solve(self, M, average, top_k = 1, inverse = False):
        if len(self.vectors) == 1:
            return [(self.sum_norm(self.vectors), self.vectors)]
        assert top_k <= len(M)*len(self.SMA_reps), """top_k must not be greater than 
        #SMA-representatives for m-VS or #vectors * #SMA-representatives for VSA"""
        candidates = []
        for SMA_rep in self.SMA_reps:
            indexed_vectors = zip(self.vectors, range(0, len(self.vectors)))
            indexed_vectors.sort(self.linear_order(SMA_rep), key = operator.itemgetter(0), reverse = not inverse)
            for vec_num in M:
                [vectors, indices] = zip(*indexed_vectors[0:vec_num])
                vec_sum_norm = self.sum_norm(vectors)
                if average:
                    vec_sum_norm /= float(len(vectors))
                candidates.append((vec_sum_norm, sorted(indices)))
        candidates.sort(key = operator.itemgetter(1))
        for i in range(0, len(candidates)-1):
            if candidates[i][1] == candidates[i+1][1]:
                candidates[i] = None
        candidates = filter(lambda x: x is not None, candidates)
        if inverse:
            return heapq.nsmallest(top_k, candidates, key = operator.itemgetter(0))
        else:
            return heapq.nlargest(top_k, candidates, key = operator.itemgetter(0))
        
    def m_VS(self, m, top_k = 1, inverse = False):
        assert m <= len(self.vectors), "m must not be greater than #vectors"
        return self.solve([m], False, top_k, inverse)
        
    def VSA(self, top_k = 1, inverse = False):
        return self.solve(range(1, len(self.vectors)+1), True, top_k, inverse)