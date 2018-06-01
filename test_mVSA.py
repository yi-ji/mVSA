from mVSA import mVSA
import unittest

class test_mVSA(unittest.TestCase):

    def setUp(self):
        self.all_vectors = [[[1, 0, 5.1], [5.2, 0, 0], [1.1, 2, 3], [2, 5, 4]], [[0.1, 2.3, 4.2, 4.6]], [[0.1, 2.3, 4.2, 4.6], [12, 3.4, 4.33, 5.555], [1, 3, 5, 6], [4.2, 5, 9, 8], [1, 4, 5.5, 6.4]], [[43, 2, 13, 3], [43, 2, 13, 3.001], [1, 2, 3, 4], [6, 1, 5, 23]]]
        self.m_VS_cases = [[([2, 3], [(10.807867504739313, [0, 3]), (10.373523991392704, [2, 3]), (9.635351576356724, [1, 3])]), ([3, 2, True], [(11.086027241532468, [0, 1, 2]), (12.918591254467339, [1, 2, 3])])], [([1, True], [(6.640783086353596, [[0.1, 2.3, 4.2, 4.6]])])], [([3, 3, True], [(24.412906422628176, [0, 2, 4]), (26.293438820359732, [0, 1, 2]), (27.14054024886019, [0, 1, 4])])], [([2, 4], [(90.13330128759293, [0, 1]), (58.395650531525035, [1, 3]), (58.39520528262573, [0, 3]), (47.508041435108645, [1, 2])])]]
        self.VSA_cases = [[([2, True], [(3.6293938887919017, [1, 2]), (3.695342413844156, [0, 1, 2])]), ([4], [(6.708203932499369, [3]), (5.403933752369657, [0, 3]), (5.2, [1]), (5.197114584074513, [0])])], [([], [(6.640783086353596, [[0.1, 2.3, 4.2, 4.6]])])], [([5], [(14.323649150967082, [1]), (13.698175060934211, [3]), (13.17599071227663, [1, 3]), (11.558406873113803, [1, 3, 4]), (11.463529125012071, [3, 4])])], [([2, True], [(5.477225575051661, [2]), (14.585952145814822, [2, 3])])]]
        self.SMA_cases = [24, 0, 120, 24]

    def tearDown(self):
        self.all_vectors = None
        self.m_VS_cases = None
        self.VSA_cases = None

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

    def test_SMA(self):
        for idx, vectors in enumerate(self.all_vectors):
            solver = mVSA(vectors)
            self.assertEqual(self.SMA_cases[idx], len(solver.SMA_reps), 'wrong SMA representatives')

if __name__ == '__main__':
    unittest.main()