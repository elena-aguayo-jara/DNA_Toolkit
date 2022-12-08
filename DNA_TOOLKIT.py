"""
Elena Aguayo Jara
Máster en Bioinformática y Biología Computacional - Universidad Autónoma de Madrid
07/12/2022
DNA Toolkit
"""
from numpy.random import choice # generate a random string
import string # generate random names with lower case letters
import random # generate random names and mutations
import os # listing files 
import time # datetime

# The class DNA_Toolkit is a class that has a constructor that sets the properties of the objects from the class
class DNA_Toolkit():
    def __init__(self):
        """Set constructor
        """
        self.chain = '' # Empty string to store the DNA sequence
        self.freqs= {} # A diccionary to store the frequencies of bases
        self.file_suffix = '' # Empty String to store the files and to be able to save them with the same name when we loading a file and modifying the sequence
         
    # Set the methods of the class DNA_Toolkit
    def __str__(self):
        """The function returns the DNA chain if it is not empty, otherwise it returns an empty string
        :return: The chain is being returned.
        """
        if self.chain == '':
            return ''
        else:
             return self.chain

    def operations_menu(self): 
        """It's a function that asks the user to select an option from a menu, and returns the option
        elected, the variable "option2". 
        """
        try:
            print(f'''
            ------ DNA TOOLKIT ------\n
                OPERATIONS MENU\n
            1.- Regenerate DNA chain
            2.- Validate DNA chain
            3.- Mutate DNA chain
            4.- Measure frequencies
            5.- Count subsequences
            6.- Syntetize the reverse and complement DNA chain
            7.- Measure %GC
            8.- mRNA transcript
            9.- Exit\n
            ------------------------
            ''')

            option2=int(input('Please, select an option: '))

        except ValueError:
            print('That was not a valid number. Please, try again...')
            option2=-1

        if option2 >=10 or option2 <=0:
            print('Invalid input. Please, try again...')
        else:
            print('The option selected is: ', option2)
            return option2

    def regenerate(self):
        """The function will re-generate the DNA chain
        :return: The DNA chain is being returned.
        """
        return self.chain
    
    def validate(self):
        """It checks if the DNA chain is valid or not. It is valid when all the nucleotides are A, C, T or G
        """
        # A list with the valid nucleotides
        nucleotides = ['A', 'C', 'G', 'T']
        # Converting the DNA sequence to uppercase.
        DNA_validate = self.chain.upper()
        
        # Validating the DNA chain.
        for nucl in DNA_validate:
            if nucl not in nucleotides:
                # If in the chain there are some characters that are not in the nucleotides list return False
                return False
        return DNA_validate 

    def mutate(self):
        """ This function mutate the DNA sequence.Thanks of the loop while, the program 
        ask you the number of mutations until you enter a number between 1 and the length of the DNA sequence
        """
        while True: 
            try:
                print('REMEMBER: The number must be between 1 and ', len(self.chain))
                num_mutations = int(input('Please, provide the number of mutations to insert in the DNA chain: '))
                
                if num_mutations <1 or num_mutations > len(self.chain):
                        print('The number is not correct. Try again...')
                else:
                    l_position = list(self.chain)
                    position = random.sample(range(0,len(self.chain)), num_mutations)

                    # Generating a random nucleotide for each position in the sequence.
                    for i in position:
                        same_nucl = l_position[i]
                        random_nucl = None
                        # Generating a random nucleotide that is different from the one in the original sequence with numpy
                        while random_nucl == None or random_nucl == same_nucl:
                            if l_position[i] == 'A':
                                random_nucl = choice(['C', 'G', 'T'])
                            elif l_position[i] == 'C':
                                random_nucl = choice(['A', 'G', 'T'])
                            elif l_position[i] == 'T':
                                random_nucl = choice(['A', 'G', 'C'])
                            elif l_position[i] == 'G':
                                random_nucl = choice(['A', 'C', 'T'])
                            else:
                                print('Mutation could not be performed.\n')
                                print('Back to the menu...')

                        l_position[i] = random_nucl
                    str_DNA_mutate=''.join(l_position)
                    # To update our DNA chain to the mutate chain and to save it in a file
                    DNA.chain = str_DNA_mutate
                    return DNA.chain  

            except ValueError:
                    print('Sorry, that was not a valid number.\n')
                    print('Back to the menu...')
            
    def measure_freqs(self):
        """It counts the number of times each base appears in the random chain.
        Uracil has been added to the list (l_bases) in order to count its frequency in the RNA.
        :return: The dictionary of the frequencies of the bases in the random chain.
        """
        l_bases = ['A', 'C', 'G', 'T', 'U']
        # Counting the number of times each base appears in the DNA sequence. 
        for d in l_bases:
            self.freqs[d] = self.chain.count(d)
        return self.freqs
      
    def count_subseq(self):
        """It's a function that asks the user to input a DNA subsequence, and it counts the
        number of times that the subsequence appears in the random chain.
        The while loop allows the user to be prompted for the sub-sequence until it is correct.
        """
        while True:
            subSeq = (input('Please, provide the DNA subsequence: '))
            resultado = self.chain.count(subSeq)
            # Converting a string to concatenate in the print
            resultado_str = str(resultado) 
            """The above code is checking if the subsequence is in the DNA chain. 
            If it is, it will return the number of times it appears. If it is not, it will return 0.
            """
            if resultado == 0:
                print('The subsequence ', subSeq, ' does not appers in the DNA chain. Try again...' )                      
            elif resultado > 1:          
                print('The subsequence ', subSeq, ' appears ', resultado_str,  ' times')
                return resultado_str

    def synt_complement(self):
        """It takes a random DNA chain and returns its reverse and complementary chain
        """   
        reverse = DNA.chain[::-1]          
        # Create an empty list to store de complementary sequence
        DNA_complement_list=[]         
        """This is a loop that iterates over the random chain, and if the base is A, it
        adds T to the DNA_complement_list, and so on. 
        If the loop find an Uracil, it will create a complementary chain of RNA
        """
        for c in DNA.chain:
            if 'U' in DNA.chain:
                if c == 'U':
                    DNA_complement_list.append('A')
                elif c == 'A':
                    DNA_complement_list.append('U')
                elif c == 'C':
                    DNA_complement_list.append('G')
                elif c == 'G':
                    DNA_complement_list.append('C')
            # If there is not an Uracil, if will create a complementary chain of DNA 
            else: 
                if c == 'A':
                    DNA_complement_list.append('T')
                elif c == 'T':
                    DNA_complement_list.append('A')
                elif c == 'C':
                    DNA_complement_list.append('G')
                elif c == 'G':
                    DNA_complement_list.append('C')

        # Converting the list in a string                    
        DNA_complement_str=''.join(DNA_complement_list)
        # To update our DNA chain to the complementary chain and to save it in a file
        self.chain = DNA_complement_str
        return reverse, DNA_complement_str
                    
    def measureGC(self):
        """It takes a sequence of DNA as input and returns the percentage of GC content in the sequence
        :param seq: the sequence of nucleotides.
        """
        GC_content = round((self.chain.count('C') + DNA.chain.count('G')) / len(self.chain) * 100)  
        return GC_content
                              
    def transcript(self):
        """It takes a DNA sequence as input and returns the corresponding RNA transcript         
        :param seq: The random DNA sequence
        """
        # Create an empty list
        l_RNA=[]
        """It is a loop that iterates over the random chain, generates the complementary DNA chain 
        changing T for U to generate the RNA transcript.
        """
        for c in self.chain:
            if c == 'A':
                l_RNA.append('U')
            elif c == 'T':
                l_RNA.append('A')
            elif c == 'C':
                l_RNA.append('G')
            elif c == 'G':
                l_RNA.append('C')
                            
        str_RNA=''.join(l_RNA)
        # To update our DNA chain to the transcript and to save it in a file
        DNA.chain = str_RNA
        return DNA.chain

# The following functions do not belong to the class, they are outside of the class
def general_menu():
    """This function prints a menu of options for the user to choose from and returns the option
    elected, the variable "option". 
    """
    option = 0
    while option < 1 or option > 7:
        try:
            print(f'''
            ------ DNA TOOLKIT ------\n
                GENERAL MENU\n
            1.- Create new DNA chain
            2.- Save DNA chain
            3.- Load DNA from disk
            4.- List all DNA chains
            5.- Delete DNA info
            6.- Operations with DNA
            7.- Exit\n
            
            -------------------------
            ''')
            # Printing the DNA chain in the menu. If the length is more than 200 bases, it will display only 200
            if len(DNA.__str__()) > 200:
                print (f"DNA chain:\n[[5' {DNA.__str__()[:200]}...]]") 
            elif DNA.__str__() == '':
                print ('[[There is not DNA chain loaded]]\n')
            else: 
                print(f"DNA chain:\n[[5' {DNA.__str__()} 3']]")

            option = int(input('\nPlease, select and option: '))
    
        except ValueError:
            print('That was not a valid number. Please, try again...')
            option=-1

        if option >=8 or option <=0:
            print('Invalid input. Please, try again...')
        else:
            print('The option selected is: ', option)
    return option
        
def random_DNA():
    """This function creates a random DNA sequence of a given length with numpy
    """
    try: 
        bases=['A','C','G','T']
        length=int(input('Please, provide the length of the DNA chain: '))   
        # Creating an empty list.
        l_DNA=[] 
        l_DNA=choice(bases,length) 
        # Concatenating the elements of the list DNA in a string.
        DNA.chain="".join(l_DNA) 
        print(f"\nThe DNA chain is:\n5' {DNA.chain} 3' ")
        return DNA.chain

    except ValueError:
        print('Sorry, that was not a valid number.\n')
        print('Back to the menu...')

def saveChain(sequence):
    """ The DNA chain that is being saved. 
    To create tha names of the files, it generates a random string of 5 characters with date and time.
    """
    # The above code is saving the DNA chain in a file if the string is empty
    if DNA.file_suffix == '':
        file_name ="".join(random.choice(string.ascii_lowercase) for i in range(5))
        datetime = time.strftime('.%d-%m-%Y.%H:%M')    
        str_datetime = str(datetime)       
        # Adding the extension .dna to the file name.
        append_suf = '.dna'
        file_suffix = (file_name + str_datetime + append_suf)
                    
        # Opening the file with the name file_suffix and it is writing the DNA chain in it.
        DNA_file = open(file_suffix, 'w')
        DNA_file.write(DNA.chain)
        DNA_file.close()
        print(f"The DNA chain 5' {DNA.__str__()} 3' is saved in your computer")
        print('Tha name of the file is', file_suffix)
        return file_suffix  
    
   # If file_suffix is not empty, it will rewrite the content to save the sequence with the same name
    else:
        DNA_file = open(DNA.file_suffix, 'w')
        DNA_file.write(DNA.__str__())
        DNA_file.close()
        print(f"The DNA chain 5' {DNA.__str__()} 3' is saved in your computer")
        print('Tha name of the file is', DNA.file_suffix)
        return DNA.file_suffix
    
def loadFile():
    """It asks the user to input a file name to load it.
    """ 
    try:  
        file_suffix= input('What file do you want to open?: ')
        f = open(file_suffix, 'r')
        # Reading the first line of the file and removing the new line character.
        DNA.chain = f.readline().rstrip('\n') 
        DNA.file_suffix = file_suffix
        print(f" The DNA sequence is:\n 5' {DNA.__str__()} 3'")
        f.close()   
        # Catching exceptions          
    except FileNotFoundError:
        print('That file does not exist.\n')
        print('Back to the menu...')
    except UnboundLocalError:
        print('The file name is not correct.\n')
        print('Back to the menu...')
    except PermissionError:
        print('The file name is not correct.\n') 
        print('Back to the menu...')
    except OSError:
        print('The file name is not correct.\n') 
        print('Back to the menu...')
    
def listFile():
    """It takes a list of files in the current directory, filters out the ones that are not .dna
    files, and prints the remaining list.
    """
    list_file = os.listdir() 
    list_dna = []
    # Checking if the file is a file and if it ends with .dna. If it does, it adds it to the list_dna
    for file in list_file:
        if os.path.isfile(os.path.join(file)) and file.endswith('.dna'):
            list_dna.append(file)
    return list_dna
       
def deleteFile():
    """It asks the user to input a file name, and if the file exists, it deletes it
    :param DNA: The DNA sequence
    """           
    try:
        delete_file= input('What file do you want to delete?: ')
        os.remove(delete_file)
    except FileNotFoundError:
        print('That file does not exist.\n')
        print('Back to the menu...')
    except UnboundLocalError:
        print('The file name is not correct.\n')
        print('Back to the menu...')
    except PermissionError:
        print('The file name is not correct.\n') 
        print('Back to the menu...')
    except OSError:
        print('That file does not exist\n')
        print('Back to the menu...')
    else:
        print('File is deleted sucessfully!')

# Below is the core of the programme. First, we define the object DNA of the class DNA_Toolkit()
DNA = DNA_Toolkit() 
exitCode=7
option=0
# A loop that allows the user to select an option from the menu while the option is not exitCode (number 7)
while option !=exitCode:  
    # The menu is defined after the while and thus the option is reevaluated in each cycle.
    option = general_menu()
    if option == 1:
        """The "if" statement is a conditional statement. It's a way of saying "if this condition is true,
        then do this". The condition is "option == 1". 
        The "==" is a comparison operator. It's a way of saying "is equal to".
        """ 
        random_DNA()     

    elif option == 2:
        # If the option is equal to 2, it is call the saveChain function
        saveChain(DNA.__str__())

    elif option == 3:
        # Calling the loadFile function
        loadFile()

    elif option == 4:
        # Calling the listFile function
        listing = listFile()
        print('LIST OF FILES: ', listing)

    elif option == 5:
        # Calling the deleteFile function
        deleteFile()
    
    elif option == 7:
        print('\nThank you. Bye!')

    elif option == 6: 
        """Checking if the DNA chain is empty. If it is, it will print a message. If it is not, it will
        set the option2 variable to 0 and the exitCode_operation to 9.
        """
        if DNA.chain == '':
            print('Please, create or load a DNA sequence first')
        else:
            option2 = 0
            exitCode_operation = 9         
            """It's a loop that allows the user to select an option from the operations_menu, 
            and returns the variable "option2" while is not exitCode_operation (number 9)
            """
            while option2 !=exitCode_operation:
                option2=DNA.operations_menu()
                    # Calling the methods of the class DNA_Toolkit
                if option2 == 1:
                    reg = DNA.regenerate()
                    print(f"The re-generate DNA chain is:\n 5' {reg} 3'")

                elif option2 == 2:
                    valid = DNA.validate()
                    # This print means that the DNA chain is valid
                    print(f"DNA sequence:\n 5' {valid} 3'\n--- Everything is OK. You can continue ---")

                elif option2 == 3:
                    no_mut = DNA.chain
                    mut = DNA.mutate()
                    print(f"\nPrevious DNA chain: 5' {no_mut} 3'\n")
                    print(f"Mutated DNA chain: 5' {mut} 3'")

                elif option2 == 4:  
                    frequences = DNA.measure_freqs()
                    print("\nUniform frequency: %d" % (len(DNA.chain)/4))
                    print(frequences)

                elif option2 == 5:
                    DNA.count_subseq()

                elif option2 == 6:
                    original = DNA.chain
                    synt = DNA.synt_complement()
                    print(f"\nThe original DNA chain is: 5' {original} 3'")
                    # Reversing the DNA chain.
                    print(f"\nThe DNA reverse is: 5' {synt[0]} 3'\n")
                    print('--- Complementary chain ---')

                    # Printing the DNA chain in a nice way when their length is <=100
                    if len(DNA.chain) <= 100:
                        print(f"\n 5' {original} 3' ")
                        print(f"    {''.join(['|' for i in range (len(DNA.chain))])} ")  
                        print(f" 3' {synt[1]} 5' ")
                    else:
                        print(f" 3' {synt[1]} 5' ")
                    
                elif option2 == 7:
                    gc = DNA.measureGC()
                    print('The GC content is:', gc, '%')

                elif option2 == 8:
                    tr = DNA.transcript()
                    print(f"RNA is  5' {tr} 3'")
                
                elif option2 == 9:
                    print('Back to the general menu...')
                




                
                
                
                






    
              
    
    
              
       

  
        
        
        
 
        
     

    



                

                
                
                
                

                
                
                

                
                        




                                





                        
