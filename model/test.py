import unittest
import glob
import subprocess

filenames = glob.glob("golden/*")
filenames = [i[i.find('/')+1:] for i in filenames]
class testFiles(unittest.TestCase):
    def test_exist(self):
        print("\n"+"=="*10+"start test_exist"+"=="*10)
        second_run = glob.glob("runs/default/output/*")
        second_run = [i[i.rfind("/")+1:] for i in second_run]
        for i in second_run:
            print("\tChecking {} ".format(i))
            self.assertTrue(i in filenames, "{} does not exist".format(i))
        print("=="*10+"test_exist compleate"+"=="*10)

    def test_compare_files(self):
        print("\n"+"=="*10+"start test_compare_files"+"=="*10)
        for comparable_file in filenames:
            golden_file_output = None
            second_file_output = None
            try:
                golden_file_data = open("golden/{}".format(comparable_file))
                golden_file_output = [line for line in golden_file_data]
                golden_file_data.close()
            except UnicodeDecodeError:
                print("\tgolden/{} is binary file".format(comparable_file))
                golden_file_data.close()
                pass # Fond non-text data
            try:
                second_file_data = open("runs/default/output/{}".format(comparable_file))
                second_file_output = [line for line in second_file_data]
                second_file_data.close()
            except UnicodeDecodeError:
                print("\truns/default/output/{} is binary file".format(comparable_file))
                second_file_data.close()
                pass # Fond non-text data

            if not second_file_output == None and not golden_file_output == None:
                index = 0
                for line in zip(second_file_output, golden_file_output):
                    self.assertEqual(line[0], line[1], "{} not the same on line {}".format(comparable_file, index))
                    index+=1
                print("\tProcessed {}".format(comparable_file))
        print("=="*10+"test_compare_files compleate"+"=="*10)


if __name__ == '__main__':
    unittest.main()
