'''
author: Elena Aguayo Jara
november 2022
Asignment: DNA Toolkit
'''
import numpy as np
def general_menu(): #generamos el menú con una función donde irán dentro el resto de operaciones
    option = 0
    while option <=1 or option>=7:
        try:
            print(f'''
                DNA TOOLKIT\n
                1.-Create new DNA chain
                2.-Save DNA chain
                3.-Load DNA from disk
                4.-List all DNA info
                5.-Delete DNA info
                6.-Operations with DNA
                7.-Exit
            ''')
            option=int(input('Please, select an option: '))
        except ValueError:
            print('That was not a valid number. Please, try again...')
            option=-1
        else:
            print('The option selected is: ', option) #FALTA PONER LA CADENA QUE HAY CARGADA
        return option
chain=0
exitCode=7
option=0
while option !=exitCode:
    option=general_menu()

#CREATE NEW DNA CHAIN
    if option == 1: 
        bases=['A','C','G','T']
        length=int(input('Please, provide the length of the DNA chain: '))
        DNA=[] #cadena creada añadida a una lista
        def random_protein(DNA):
            for i in bases: #bucle for para que recorra cada nucleótido y cree la cadena
                DNA=np.random.choice(bases,length) #crea cadena aleatoria
                chain="".join(DNA) #concatenar
           
            print(chain)
            return chain
        
        random_protein(DNA)

#FILE OPERATIONS WITH DNA CHAIN
#SAVE
    def save(chain):
        if option == 2:
            file_DNA= open('DNA_chain.dna', 'w')
            for i in range(len(chain)):
                file_DNA.write(chain[i]+',')
            file_DNA.close()
    save(chain)

general_menu()