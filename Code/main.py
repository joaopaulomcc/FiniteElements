from job import *

print("######################################################################")
print("##                    SOLVER - Joao Paulo                           ##")
print("######################################################################")

# Main loop

while True:

    job()
    quit_yn = input("Do you want to quit the program? (Y/N) ").upper()

    if quit_yn == 'Y':
        break
