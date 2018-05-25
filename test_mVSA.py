from mVSA import mVSA
import unittest

class test_mVSA(unittest.TestCase):

    def setUp(self):
        self.all_vectors = [[[1,0,5.1], [5.2,0,0], [1.1,2,3], [2,5,4]], [[0.1, 2.3, 4.2, 4.6]]]
        self.m_VS_cases = [[([2, 3], [(10.807867504739313, [0, 3]), (10.373523991392704, [2, 3]), (9.635351576356724, [1, 3])]), ([3, 2, True], [(11.086027241532468, [0, 1, 2]), (12.918591254467339, [1, 2, 3])])], [([1, True], [(6.640783086353596, [[0.1, 2.3, 4.2, 4.6]])])]]
        self.VSA_cases = [[([2, True] ,[(3.6293938887919017, [1, 2]), (3.695342413844156, [0, 1, 2])]), ([4], [(6.708203932499369, [3]), (5.403933752369657, [0, 3]), (5.2, [1]), (5.197114584074513, [0])])], [([], [(6.640783086353596, [[0.1, 2.3, 4.2, 4.6]])])]]

    def tearDown(self):
        self.cases = None

    def test_m_VS(self):
        for idx, vectors in enumerate(self.all_vectors):
            solver = mVSA(vectors)
            for case in self.m_VS_cases[idx]:
                answer = solver.m_VS(*(case[0]))
                truth = case[1]
                self.assertEqual(len(answer), len(truth), 'invalid answer')
                for i in range(0, len(truth)):
                    self.assertEqual(answer[i], truth[i], 'wrong answer')

    def test_VSA(self):
        for idx, vectors in enumerate(self.all_vectors):
            solver = mVSA(vectors)
            for case in self.VSA_cases[idx]:
                answer = solver.VSA(*(case[0]))
                truth = case[1]
                self.assertEqual(len(answer), len(truth), 'invalid answer')
                for i in range(0, len(truth)):
                    self.assertEqual(answer[i], truth[i], 'wrong answer')

if __name__ == '__main__':
    unittest.main()