import unittest
import glob
import subprocess

# Grabs the files in golden
filenames = glob.glob("golden/*")
# cleans the filenames
filenames = [i[i.find('/')+1:] for i in filenames]

class testFiles(unittest.TestCase):

    # checks if all the files in golden are also in the default output
    def test_exist(self):
        print("\n"+"=="*10+"start test_exist"+"=="*10)

        # grabs the files in default/output
        second_run = glob.glob("runs/default/output/*")
        # cleans the file names
        second_run = [i[i.rfind("/")+1:] for i in second_run]

        # checks that they all exist
        for i in second_run:
            print("\tChecking {} ".format(i))
            self.assertTrue(i in filenames, "{} does not exist".format(i))
        # print compleate
        print("=="*10+"test_exist compleate"+"=="*10)

    # checks if the files have the same values
    def test_compare_files(self):
        print("\n"+"=="*10+"start test_compare_files"+"=="*10)
        all_error = ""
        # for every file in both the golden and default runs
        for comparable_file in filenames:
            error = ""
            golden_file_output = None
            second_file_output = None

            # grab the file if they exist and save all the values
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

            # If both the files have been found and had lines in them
            if not second_file_output == None and not golden_file_output == None:
                # check if everyline is good
                line_num = 0
                for lines in zip(second_file_output, golden_file_output):
                    second_file_line = lines[0].split()
                    golden_file_line = lines[1].split()
                    for words in zip(second_file_line, golden_file_line):
                        if words[0][0].isdigit() and words[1][0].isdigit():
                            number1 = float(words[0])
                            number2 = float(words[1])
                            if abs(number1 - number2) > (1/10**8):
                                error = error + "\nline num: {:<20} output:{:<20} Expected:{}".format(line_num,number1, number2)
                        else:
                            if words[0] != words[1]:
                                    error = error + "\nline num: {:<20} output:{:<20} Expected:{}".format(line_num,words[0][:-2], words[1][:-2])
                    line_num= line_num+1
            if error!="":
                error = "\n\n"+ comparable_file +"\n" + error
                all_error = error + all_error

        print(all_error)
        self.assertEqual(all_error,"")
        print("\tProcessed {}".format(comparable_file))
        print("=="*10+"test_compare_files compleate"+"=="*10)


if __name__ == '__main__':
    unittest.main()
