import math
import itertools

class Seq_alignment(): 
    '''
    This class contains a total of 7 methods and an initializing constructor.
    Let us see what each function is doing in this program:
        1) __init__ - This initializes the variables seq_1,seq_2,l1,l2 through the class object "construct_sa".
        2) permutation_tuple_duplicates - This creates two two-dimensional arrays of the different permutations of the respective nucleic acid strings (including the gaps).
        3) remove_duplicates - This function removes the duplicates in the two arrays created in the previous function.
        4) twod_to_oned - This function converts the two dimensional arrays into one dimensional arrays for further processing.
        5) sequence_array - This function removes the "_" character from each string of the two one dimensional arrays to compare the resultant strings with the two respective original sequences for futher elimination.
        6) deletion - This function removes the strings from the two one-dimensional arrays whose sequences (without the gaps: "_"), don't match the original two sequences. 
        7) oned_to_twod - This function creates a new one-dimensional array which contains all possible combinations of strings from both the first and second string. Next, it converts this one-dimensional array into a two-dimensional array with each sub-array of length 2 to compare each permutation of the two sequences (with/without gaps).
        8) underscore_comp - This function compares the two strings in each sub-array of the two-dimensional array character by character to check if a gap doesn't align with another gap because it doesn't make sense.
    '''
    
    def __init__(self,seq_1,seq_2,l1,l2):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.l1 = l1
        self.l2 = l2
    
    def permutation_tuple_duplicates ():
        string_1=''
        string_2=''
        string_1_array_dup=[]
        string_2_array_dup=[]
        for i in range(l2-l1,l2+1):
            if l1<l2:
                string_1 = seq_1+i*"_"
                string_2 = seq_2+(i-1)*"_"
                string_1_array_dup.append(["".join(perm) for perm in itertools.permutations(string_1)])
                string_2_array_dup.append(["".join(perm) for perm in itertools.permutations(string_2)])
            elif l1==l2:
                string_1 = seq_1+i*"_"
                string_2 = seq_2+i*"_"
                string_1_array_dup.append(["".join(perm) for perm in itertools.permutations(string_1)])
                string_2_array_dup.append(["".join(perm) for perm in itertools.permutations(string_2)])
        return (string_1_array_dup,string_2_array_dup)
    
    def remove_duplicates(dup_1, dup_2):
        string_1_array_temp=[]
        string_2_array_temp=[]
        for i in range(0,l1+1):
                string_1_dup = [*set(dup_1[i])]
                string_2_dup = [*set(dup_2[i])]
                string_1_array_temp.append(string_1_dup)
                string_2_array_temp.append(string_2_dup)
        return (string_1_array_temp,string_2_array_temp)
    
    def twod_to_oned (non_dup_1_2d,non_dup_2_2d):
        non_dup_1_1d = []
        non_dup_2_1d = []
        for i in non_dup_1_2d:
            for j in i:
                non_dup_1_1d.append(j)
        for k in non_dup_2_2d:
            for l in k:
                non_dup_2_1d.append(l)
        return (non_dup_1_1d,non_dup_2_1d)
    
    def sequence_array (seq_arr_1,seq_arr_2):
        alpha_array_1 = []
        alpha_array_2 = []
        for i in seq_arr_1:
            str_1 = str(i)
            sequence_1 = str_1.replace("_","")
            alpha_array_1.append(sequence_1)
       
        for j in seq_arr_2:
            str_2 = str(j)
            sequence_2 = str_2.replace("_","")
            alpha_array_2.append(sequence_2)
        
        return (alpha_array_1,alpha_array_2)
    
    def deletion (seq_array_1,seq_array_2,string_array_1,string_array_2):
        del_array_1=[]
        del_array_2=[]
        len_array_1 = len(seq_array_1)
        len_array_2 = len(seq_array_2)
        for i in range(0,len_array_1):
            if seq_array_1[i] != seq_1:
                del_array_1.append(string_array_1[i])
            else:
                continue
        for j in range(0,len_array_2):
            if seq_array_2[j] != seq_2:
                del_array_2.append(string_array_2[j])
            else:
                continue
        final_arr_1 = []
        final_arr_2 = []
        final_arr_1 = list(set(string_array_1)-set(del_array_1))
        final_arr_2 = list(set(string_array_2)-set(del_array_2))
        
        return (final_arr_1,final_arr_2)
    
    def oned_to_twod(final_seq_one,final_seq_two):
        oned_arr_seq = []
        len_1 = len(final_seq_one)
        len_2 = len(final_seq_two)
        for i in range (0,len_1):
            for j  in range(0,len_2):
                oned_arr_seq.extend((final_seq_one[i],final_seq_two[j]))
        
        twod_arr_seq = [oned_arr_seq[k:k + 2] for k in range(0, len(oned_arr_seq), 2)]
        
        return twod_arr_seq
       
    def underscore_comp(final_seq_arr):
        len_final = len(final_seq_arr)
        final_seq_arr_del = []
        for i in range(0,len_final):
            a = final_seq_arr[i][0]
            b = final_seq_arr[i][1]
            for x,y in zip(a,b):
                if x==y and x=="_":
                    final_seq_arr_del.append(final_seq_arr[i])
                    break
                else:
                    continue
        
        for j in final_seq_arr_del:
            final_seq_arr.remove(j)
        
        final_del_arr = []
        
        for i in range(0,len(final_seq_arr)):
            if len(final_seq_arr[i][0]) != len(final_seq_arr[i][1]):
                final_del_arr.append(final_seq_arr[i])
            else:
                continue
        
        for j in final_del_arr:
            final_seq_arr.remove(j)
         
        return (final_seq_arr,len(final_seq_arr))    
                
print("Please enter two sequences of nucleic acid (either DNA or RNA): ") #Assume l2 is greater than l1
seq_1 = input()
seq_2 = input()
l1 = len(seq_1)
l2 = len(seq_2)

construct_sa = Seq_alignment(seq_1, seq_2, l1, l2)

perm_tuple = Seq_alignment.permutation_tuple_duplicates()

string_1_array_duplicate = perm_tuple[0]
string_2_array_duplicate = perm_tuple[1]

non_duplicates = Seq_alignment.remove_duplicates(string_1_array_duplicate,string_2_array_duplicate)

string_1_array_non_duplicate = non_duplicates[0]
string_2_array_non_duplicate = non_duplicates[1]

non_dup_1d_tup = Seq_alignment.twod_to_oned(string_1_array_non_duplicate,string_2_array_non_duplicate)

string_1_array_non_duplicate_1d = non_dup_1d_tup[0]
string_2_array_non_duplicate_1d = non_dup_1d_tup[1]

sequence_tuple = Seq_alignment.sequence_array(string_1_array_non_duplicate_1d,string_2_array_non_duplicate_1d)
string_1_array_sequences = sequence_tuple[0]
string_2_array_sequences = sequence_tuple[1]

string_array_temp_tuple = Seq_alignment.deletion(string_1_array_sequences,string_2_array_sequences,string_1_array_non_duplicate_1d,string_2_array_non_duplicate_1d)
final_seq_1 = string_array_temp_tuple[0]
final_seq_2 = string_array_temp_tuple[1]

string_array_temp_2d_value = Seq_alignment.oned_to_twod(final_seq_1,final_seq_2)

final_arr_and_length = Seq_alignment.underscore_comp(string_array_temp_2d_value)
print("These are all the possbile alignments: ")
print(final_arr_and_length[0])
print("\n\n")
print("The total number of alignments possible for these sequences are :",final_arr_and_length[1])
