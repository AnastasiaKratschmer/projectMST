import random as random
import Bio
import numpy as np
import sys
import os
import networkx as nx
import random as random
from tqdm import tqdm # loading bar
from utils_copy import linear_C, get_cost_2, get_sequence_string, parse_fasta_multiple, create_score_matrix, write_alignments_to_file, linear_backtrack, fill_graph,new_sp_approxi_combi
from utils_copy import convert_to_desired_format_nr_version, compute_cost, my_traversal_simply, extend_alignment_chaos, find_min_span_edges_testing, parse_fasta_multiple_remove_n
import timeit
from utils_copy import al_integrity_testt
from old_for_testing.sp_approx import sp_approx
from old_for_testing.utils import *
from utils_copy import sum_of_column


score_matrix={'a': {'a': 0, 'c': 5, 'g': 2, 't': 5}, 'c': {'a': 5, 'c': 0, 'g': 5, 't': 2}, 'g': {'a': 2, 'c': 5, 'g': 0, 't': 5}, 't': {'a': 5, 'c': 2, 'g': 5, 't': 0}}
gap_cost=5

def find_highest_second_element(d, first_element): #to find the highest col nr already in the MSA
    highest_second_element = 0
    for key, value in d.items():
        if value[0] == first_element:
            second_element = value[1]
            if second_element > highest_second_element:
                highest_second_element = second_element
    return highest_second_element

def finish_running_through_MSA1(MSA1,MSA2,j,i,l):
    how_many_cols_in_MSA2=len(MSA2[0])
    gaps_to_add=['-']*how_many_cols_in_MSA2
    new_col=list(MSA1[j])+gaps_to_add
    return new_col

def finish_running_through_MSA2(MSA1,MSA2,j,i,l):
    how_many_cols_in_MSA1=len(MSA1[0])
    gaps_to_prepend=['-']*how_many_cols_in_MSA1
    new_col=gaps_to_prepend+list(MSA2[l])
    return new_col

def equality_here_test(a,b):
    equality=False
    if a=='-' and b=='-' or a!='-' and b!='-':
        equality=True
    return equality

def new_identify_merge_case_compact(guide_PA,MSA1,MSA2,node1,node2,MSA_position_overview_dict,i,j,l):
    node1_position=MSA_position_overview_dict[node1]
    node1_col=node1_position[1]
    node2_position=MSA_position_overview_dict[node2]
    node2_col=node2_position[1]
    case=None
    move_on_in_guide_flag=False
    a=guide_PA[i][0]
    b=guide_PA[i][1]
    c= MSA1[j][node1_col]
    d=MSA2[l][node2_col]
    look_forward_flag=False

    if equality_here_test(a,c)==True and equality_here_test(b,d)==True: case=3 #a==c and b==d
    elif equality_here_test(a,c)==False and equality_here_test(b,d)==False:
        if equality_here_test(b,c)==True:
            if MSA1[j+1][node1_col]==guide_PA[i][0]:
                look_forward_flag="col1"
            if MSA2[l+1][node2_col]==guide_PA[i][1]:
                look_forward_flag="col2"
        case=4
    elif equality_here_test(a,c)==True:#a==c
        if d!='-':
            case= 2 #was orginially 1 2 omg just seeing what will happen if I put it to 1. #truing to reverse to see
            move_on_in_guide_flag=True
        else: case= 1
    elif equality_here_test(a,c)==False: #a!=c
        if a=='-':
            case= 1 #was originally 1 omg!, just checking what will happen if I put it to 2 #reversing to see what will happen...
            move_on_in_guide_flag=True
        else: case=2
    else: sys.exit()
    guide_was=guide_PA[i]
    MSA2_was=MSA2[l]
    return case,move_on_in_guide_flag,guide_was,MSA2_was, look_forward_flag

def alt_handle_case1(MSA1,MSA2,i,j,l):
    how_many_cols_in_MSA1=len(MSA1[0])
    gaps_to_prepend=['-']*how_many_cols_in_MSA1
    new_col=gaps_to_prepend+list(MSA2[l])
    return new_col #after this, increase only l
def alt_handle_case2(MSA1,MSA2,i,j,l):
    how_many_cols_in_MSA2=len(MSA2[0])
    gaps_to_append=['-']*how_many_cols_in_MSA2
    col_now=MSA1[j]
    col_new=list(MSA1[j])+gaps_to_append
    return col_new #after this, increase only j
def alt_handle_case3(MSA1,MSA2,i,j,l):
    col_new=list(MSA1[j])+list(MSA2[l])
    return col_new #after this, increase i, j and l
def alt_handle_case4(MSA1,MSA2,i,j,l):
    #col_new=list(MSA1[j])+list(MSA2[l])
    how_many_cols_in_MSA1=len(MSA1[0])
    how_many_cols_in_MSA2=len(MSA2[0])
    #make first column!
    first_col_gaps=['-']*how_many_cols_in_MSA2
    first_col=list(MSA1[j])+first_col_gaps
    #make second colum!
    second_col_gaps=['-']*how_many_cols_in_MSA1
    second_col=second_col_gaps+list(MSA2[l])
    return first_col,second_col #increade l and j after 


def alt_alt_merge_united(guide_PA,MSA_list,in_which_MSA_is_it,node1,node2):
    print(type(in_which_MSA_is_it))
    print(in_which_MSA_is_it)
    print(type(node1))
    print(type(node2))
    MSA_of_node1=in_which_MSA_is_it[node1][0]
    MSA_of_node2=in_which_MSA_is_it[node2][0]
    spot_of_node1=in_which_MSA_is_it[node1][1]
    spot_of_node2=in_which_MSA_is_it[node2][1]
    
    print("the node1 is right now in the alignment nr "+str(MSA_of_node1))
    print("the node2 is right now in the alignment nr "+str(MSA_of_node2))
    MSA1=MSA_list[int(MSA_of_node1)]
    MSA2=MSA_list[int(MSA_of_node2)]
    j=0 
    i=0
    l=0
    MSA_new=[]
    while i<=(len(guide_PA)-1) and j<=(len(MSA1)-1) and l<=(len(MSA2)-1):
        new_col=[]
        print("I don't have the case identified yet, but MSA1 is: " + str(MSA1[j]) + " MSA2 is: " + str(MSA2[l]) + " ,and the guide is: " + str(guide_PA[i]))
        print("nodes pointing to cols: " + str(spot_of_node1) + " , " + str(spot_of_node2))
        print("i,j,l: "+str(i)+','+str(j)+','+str(l))
        print("max for those should be: "+str(len(guide_PA))+','+str(len(MSA1))+','+str(len(MSA2)))
        case,move_on_in_guide_flag,guide_was,MSA2_was,look_forward_flag=new_identify_merge_case_compact(guide_PA,MSA1,MSA2,node1,node2, in_which_MSA_is_it,i,j,l)
        print("flag: "+str(move_on_in_guide_flag))
        print("case:"+str(case))
        print("guide_was: "+str(guide_was))
        if case==1:
            col=alt_handle_case1(MSA1,MSA2,i,j,l)
            MSA_new.append(col) 
            if move_on_in_guide_flag==True:
                print("the damn flag was true!")
                i+=1
            l+=1
            print("I just made the new col: "+ str(col))
            print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==2:
            col=alt_handle_case2(MSA1,MSA2,i,j,l)
            MSA_new.append(col)
            if move_on_in_guide_flag==True:
                print("the flag was true dammit")
                i+=1
            j+=1
            print("I just made the new col: "+ str(col))
            print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==3:
            col=alt_handle_case3(MSA1,MSA2,i,j,l)
            MSA_new.append(col)
            i+=1
            j+=1
            l+=1
            print("I just made the new col: "+ str(col))
            print("MSA_new now has len: "+ str(len(MSA_new)))
            new_col.append(col)
        if case==4:
            col1,col2=alt_handle_case4(MSA1,MSA2,i,j,l)
            if look_forward_flag=='col1':
                MSA_new.append(col1)
                new_col.append(col1)
                j+=1
            elif look_forward_flag=='col2':
                MSA_new.append(col2)
                new_col.append(col2)
                l+=1
            else:
                MSA_new.append(col1)
                MSA_new.append(col2)
                j+=1
                l+=1
                print("I just made the new col: "+ str(col1))
                print("I just made the new col: "+ str(col2))
                print("MSA_new now has len: "+ str(len(MSA_new)))
                new_col.append(col1)
                new_col.append(col2)
        for object in new_col:
            print(object)
            second_element_in_merged_pair=object[(-(len(MSA2_was))+spot_of_node2)]
            selected_stuff=[object[int(spot_of_node1)],second_element_in_merged_pair]
            print("selected stuff from col made: ")
            print(selected_stuff)
            if case==1 and move_on_in_guide_flag==False or case==2 and move_on_in_guide_flag==False or case==4:
                if guide_was==selected_stuff: print("ALERT!")
            if guide_was!=selected_stuff and selected_stuff!=['-', '-']: print("ALARM!")
            if case==1 or case==2 and move_on_in_guide_flag==True:
                if guide_was!=selected_stuff and selected_stuff!=['-', '-']: print("WARNING! i is gonna be increased, maybe unrightfully")
            
    if j <= (len(MSA1)-1):
        while j<= (len(MSA1)-1):
            column=finish_running_through_MSA1(MSA1,MSA2,j,i,l)
            j+=1
            #l+=1
            MSA_new.append(column)
            print("I just made the new col: "+ str(column))
            print("MSA_new now has len: "+ str(len(MSA_new)))
    if l<= (len(MSA2)-1):
        while l<= (len(MSA2)-1):
            column=finish_running_through_MSA2(MSA1,MSA2,j,i,l)
            l+=1
            MSA_new.append(column)
            print("I just made the new col: "+ str(column))
            print("MSA_new now has len: "+ str(len(MSA_new)))
    if i<=(len(guide_PA)-1):
        print('yikes, the strings have run out but the guide has not')
        i+=1
    print("\n\n")
    return MSA_new