import sys
import unittest
#print("I search: \n")
#print(sys.path)
from all_functions_for_testing_mid_aug import *

class TestMSTfuncs(unittest.TestCase):
    def test_convert_to_desired_format_number_version(self): #does the convestion work on short cases?
        expected = [['', '1', '0', '1'], ['', '2', '0', '2'], ['', '3', '1', '2']]
        input_matrix = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        result = convert_to_desired_format_nr_version(input_matrix)
            
        # Use numpy's array_equal function to compare arrays
        self.assertTrue(np.array_equal(result, expected))

    def test_find_min_span_edges_testing_short(self): #does the finding of min span edges work on short cases where the edges are already sorted by lenght?
        expected=[['*', '1', '0', '1'],['*', '2', '0', '2'],['','3','0','2'],['', '3', '1', '2']]
        input_matrix=np.array([['', '1', '0', '1'], ['', '2', '0', '2'],['','3','0','2'], ['', '3', '1', '2']])
        result=find_min_span_edges_testing(input_matrix)
        error_message = f"Expected:\n{expected}\n\nActual:\n{result}"
        self.assertTrue(np.array_equal(result, expected), error_message)

    def test_find_min_span_edges_testing_longer_sort_yourself(self):#does it also work if it needs to sort cases by lenght itself?
        expected=[['*','0','0','1'],['', '1', '0', '1'], ['*', '2', '0', '2'],['','3','0','2'], ['', '3', '1', '2'],['*','5','2','3'],['','7','3','0'],['*','14','4','3']]
        input_matrix=np.array([['', '1', '0', '1'], ['', '2', '0', '2'],['','3','0','2'], ['', '3', '1', '2'],['','5','2','3'],['','7','3','0'],['','14','4','3'],['','0','0','1']])
        result=find_min_span_edges_testing(input_matrix)
        error_message = f"Expected:\n{expected}\n\nActual:\n{result}"
        self.assertTrue(np.array_equal(result, expected), error_message)

   # def test_find_min_span_edges_testing_longer_sort_yourself(self):#does it also work if it need to sort cases by lenght itself,a and sort node names? #no but does it ever need sorting?
    #    expected=[['*','0','0','1'],['', '1', '0', '1'], ['*', '2', '0', '2'],['','3','0','2'], ['', '3', '1', '2'],['*','5','2','3'],['','7','3','0'],['*','14','4','3']]
     #   input_matrix=np.array([['', '1', '0', '1'], ['', '2', '0', '2'],['','3','0','2'], ['', '3', '1', '2'],['','5','2','3'],['','7','3','0'],['','14','4','3'],['','0','1','0']])
      #  result=find_min_span_edges_testing(input_matrix)
       # error_message = f"Expected:\n{expected}\n\nActual:\n{result}"
        #self.assertTrue(np.array_equal(result, expected), error_message)
        
    def test_get_visiting_order_for_graph(self):#this was thought as the way to get the visiting order, but now I just use it for graph structure!!
        input1=np.array([['*','0','0','1'],['', '1', '0', '1'],['*', '2', '0', '2'],['','3','0','2'],['', '3', '1', '2'],['*','5','2','3'],['','7','3','0'],['*','14','4','3']])
        input2='0'
        graph=get_visiting_order(input1,input2)[1]
        self.assertEqual(graph.number_of_nodes(), 5)
        self.assertEqual(graph.number_of_edges(), 4)


if __name__ == '__main__':
    unittest.main()
