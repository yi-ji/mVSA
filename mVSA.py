from ppl import *
from itertools import imap, groupby
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
        
    def is_close(self, a, b):
        min_abs = min(abs(a), abs(b))
        return min_abs == 0 or abs(a-b)/min_abs < 1e-7
        
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
    
    def support_point(self, generators):
        for generator in generators:
            if generator.is_point():
                return generator.coefficients()
        return None
    
    def intersect(self, generators, vec):
        for generator in generators:
            if generator.is_ray():
                ray_coord = generator.coefficients()
                if self.is_close(sum([ray_coord[i]*vec[i] for i in range(0, self.dim)]), 0):
                    return True
        return False
    
    def SMA_representative(self, vec, vectors, pos_ineq, neg_ineq):
        cs = Constraint_System()
        if pos_ineq:
            cs.insert(self.hyperplane(vec) > 0)
        else:
            cs.insert(self.hyperplane(vec) < 0)
        for v in vectors:
            if vec == v:
                continue
            v_hp = self.hyperplane(v)
            cs_pos = Constraint_System(cs)
            cs_neg = Constraint_System(cs)
            cs_pos.insert(v_hp > 0)
            cs_neg.insert(v_hp < 0)
            generators_pos = NNC_Polyhedron(cs_pos).minimized_generators()
            generators_neg = NNC_Polyhedron(cs_neg).minimized_generators()
            if self.intersect(generators_pos, vec):
                if self.intersect(generators_neg, vec):
                    if neg_ineq:
                        cs.insert(v_hp > 0)
                    else:
                        cs.insert(v_hp < 0)
                else:
                    cs.insert(v_hp > 0)
            else:
                cs.insert(v_hp < 0)
        generators = NNC_Polyhedron(cs).minimized_generators()
        return self.support_point(generators)

    def SMA_representatives(self, vectors):
        SMA_reps = []
        for vec in vectors:
            SMA_reps += [self.SMA_representative(vec, vectors, True, True), 
                         self.SMA_representative(vec, vectors, True, False),
                         self.SMA_representative(vec, vectors, False, True),
                         self.SMA_representative(vec, vectors, False, False)]
        SMA_reps = list(imap(operator.itemgetter(0), groupby(sorted(SMA_reps)))) # make unique
        return SMA_reps
        
    def m_VS(self, m, top_k = 1, inverse = False):
        assert m <= len(self.vectors), "m must not be greater than #vectors"
        if len(self.vectors) == 1:
        	return [(self.sum_norm(self.vectors), self.vectors)]
        assert top_k <= len(self.SMA_reps), "top_k must not be greater than #SMA-representatives"
        candidates = []
        for SMA_rep in self.SMA_reps:
            def linear_order(vec1, vec2):
                v = sum([(vec1[i]-vec2[i])*SMA_rep[i] for i in range(0, self.dim)])
                return 1 if v > 0 else -1 if v < 0 else 0
            indexed_vectors = zip(self.vectors, range(0, len(self.vectors)))
            indexed_vectors.sort(linear_order, key = operator.itemgetter(0), reverse = not inverse)
            [vectors, indices] = zip(*indexed_vectors[0:m])
            vec_sum_norm = self.sum_norm(vectors)
            candidates.append((vec_sum_norm, sorted(indices)))
        candidates = list(imap(operator.itemgetter(0), groupby(sorted(candidates))))
        if inverse:
            return heapq.nsmallest(top_k, candidates, key = operator.itemgetter(0))
        else:
            return heapq.nlargest(top_k, candidates, key = operator.itemgetter(0))
        
    def VSA(self, top_k = 1, inverse = False):
    	if len(self.vectors) == 1:
        	return [(self.sum_norm(self.vectors), self.vectors)]
        assert top_k <= len(self.vectors)*len(self.SMA_reps), "top_k must not be greater than #vectors * #SMA-representatives"
        candidates = []
        for SMA_rep in self.SMA_reps:
            def linear_order(vec1, vec2):
                v = sum([(vec1[i]-vec2[i])*SMA_rep[i] for i in range(0, self.dim)])
                return 1 if v > 0 else -1 if v < 0 else 0
            indexed_vectors = zip(self.vectors, range(0, len(self.vectors)))
            indexed_vectors.sort(linear_order, key = operator.itemgetter(0), reverse = not inverse)
            for vec_num in range(1, len(self.vectors)):
                [vectors, indices] = zip(*indexed_vectors[0:vec_num])
                vec_sum_norm_avg = self.sum_norm(vectors)/float(len(vectors))
                candidates.append((vec_sum_norm_avg, sorted(indices)))
        candidates = list(imap(operator.itemgetter(0), groupby(sorted(candidates))))
        if inverse:
            return heapq.nsmallest(top_k, candidates, key = operator.itemgetter(0))
        else:
            return heapq.nlargest(top_k, candidates, key = operator.itemgetter(0))
