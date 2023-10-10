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


def make_MSA_position_overview_dict(res_matrix):#initializes the overview dictionary
    in_which_MSA_is_it={}
    names= np.unique(res_matrix[:, 2:])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    print(in_which_MSA_is_it)


def find_highest_second_element(d, first_element): #to find the highest col nr already in the MSA
    highest_second_element = 0
    for key, value in d.items():
        if value[0] == first_element:
            second_element = value[1]
            if second_element > highest_second_element:
                highest_second_element = second_element
    return highest_second_element

def update_MSA_position_overview_dict(MSA_position_overview_dict,node1,node2): #may be somehow inadequate!
    print(node1,node2)
    node1_origin_MSA=MSA_position_overview_dict[node1][0]
    node2_origin_MSA=MSA_position_overview_dict[node2][0]
    node1_origin_MSA_col=MSA_position_overview_dict[node1][1]
    node2_origin_col=MSA_position_overview_dict[node2][1]
    #now I have the original positions in from the position overview dict!
    nr_we_got_to=find_highest_second_element(MSA_position_overview_dict,MSA_position_overview_dict[node1])
    MSA_position_overview_dict[node2]=[node1_origin_MSA,(nr_we_got_to+1)] #update the pos_dict to show that node2 changed MSA and make the col one bigger than the max
    return MSA_position_overview_dict

def get_MSA_pos(MSA_position_overview_dict, node):
    value=MSA_position_overview_dict[node]
    MSA_pos=value[0]
    MSA_pos_col=value[1]
    return MSA_pos,MSA_pos_col

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

def alt_alt_merge_united(guide_PA,MSA_list,in_which_MSA_is_it,node1,node2):
    print(type(in_which_MSA_is_it))
    print(in_which_MSA_is_it)
    print(type(node1))
    print(type(node2))
    MSA_of_node1=in_which_MSA_is_it[node1][0]
    MSA_of_node2=in_which_MSA_is_it[node2][0]
    spot_of_node1=in_which_MSA_is_it[node1][1]
    spot_of_node2=in_which_MSA_is_it[node2][1]
    print(MSA_of_node1)
    print(MSA_of_node2)
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

def new_assembly(seqs,score_matrix,gapcost):
    # Make a matrix to hold pairwise alignment costs for all alignment combinations!
    matrix = np.full((len(seqs), len(seqs)), np.nan)
    # Loop over all pairs
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            matrix[i, j] = get_cost_2(linear_C(gap_cost, score_matrix, seq1, seq2))
        print("Here comes the distance matrix produced by the alignments: \n")
        print(matrix)
    matrix_for_MST=matrix #copy the matrix, so that we can keep the old matrix and make a changed version to the "pseudomatrix" version
    matrix_for_MST=convert_to_desired_format_nr_version(matrix_for_MST) #making the "pseudomatrix"
    print("matrix for MST: "+str(matrix_for_MST))
    min_span_edges_res=find_min_span_edges_testing(matrix_for_MST)
    print("min span edges: "+str(min_span_edges_res))
    in_which_MSA_is_it={}
    names= np.unique(min_span_edges_res[:, 2:])
    in_which_MSA_is_it ={name: [int(name),0] for name in names}
    MSA_list=[[[char] for char in seq] for seq in seqs]
    who_aligned_to_who=[]
    k=0
    for row in min_span_edges_res:
        k+=1
        if row[0]=="*":
            node1=row[2]
            node2=row[3]
            print("\n \n \n these are the nodes for the iteration "+ str(k))
            print(node1,node2)
            print("which correspond to these strings I align: "+ str(seqs[int(node1)])+" , "+str(seqs[int(node1)]) )
            who_aligned_to_who.append([node1,node2])
            cost=linear_C(gap_cost,score_matrix,seqs[int(node1)],seqs[int(node2)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(node1)], seqs[int(node2)], cost, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            A = [list(e) for e in zip(alignment1,alignment2)]
            print("original alignment, which is gonna be the guide"+str(A))
            united_MSA_new=alt_alt_merge_united(A,MSA_list,in_which_MSA_is_it,node1,node2)
            print("here we have the union: "+str(united_MSA_new))
            which_spot_in_MSA_list_to_update=min(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
            which_spot_in_MSA_list_to_remove=max(in_which_MSA_is_it[node1][0],in_which_MSA_is_it[node2][0])
            MSA_list[which_spot_in_MSA_list_to_update]=united_MSA_new
            MSA_list.pop(which_spot_in_MSA_list_to_remove)
            companions_to_update=[]
            for key, value in in_which_MSA_is_it.items():
                how_many_cols_already_in_MS1=find_highest_second_element(in_which_MSA_is_it,which_spot_in_MSA_list_to_update)
                if value[0]==which_spot_in_MSA_list_to_remove:
                    companions_to_update.append(key)
                if value[0]>which_spot_in_MSA_list_to_remove:
                    value[0]=(value[0]-1)
            for companion in companions_to_update:
                print(in_which_MSA_is_it[companion])
                in_which_MSA_is_it[companion][0]=which_spot_in_MSA_list_to_update
                col_of_element_in_old_MSA2=in_which_MSA_is_it[companion][1]
                in_which_MSA_is_it[companion][1]=(how_many_cols_already_in_MS1+col_of_element_in_old_MSA2+1)
        print("this is the updated dict after it "+str(k)+": "+str(in_which_MSA_is_it))
    print(MSA_list)
    #integrity check part 1, to check if each string is the same before and after, except for gaps.
    for i,seq in enumerate(seqs):
        col_to_extract=in_which_MSA_is_it[str(i)][1]
        j=0
        new_str_with_gaps=[]
        new_str_no_gaps=[]
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][col_to_extract]
            j+=1
            new_str_with_gaps.append(found)
            if found !='-':
                new_str_no_gaps.append(found)
        new_str_no_gaps=''.join(new_str_no_gaps)
        new_str_with_gaps=''.join(new_str_with_gaps)
        if new_str_no_gaps==seq:
            print("integrity check 1 passed for seq "+str(i))
        else:
            print("Yikes, integrity check 1 did not pas for seq "+str(i)+". constrast( new, orig): \n"+str(new_str_no_gaps)+"\n"+str(seq))
     #part 2 lol, are the alignments preserved, expect for gaps???
    print("structure of who_aligned_to_who: ")
    print(who_aligned_to_who)
    for element in who_aligned_to_who:
        seq1_nr=element[0]
        seq2_nr=element[1]
        print("seq1_nr and seq2_nr are: "+str(seq1_nr)+" , "+str(seq2_nr))
        pos_in_MSA_seq1=in_which_MSA_is_it[seq1_nr][1]
        pos_in_MSA_seq2=in_which_MSA_is_it[seq2_nr][1]
        print("pos_in_MSA_seq1,pos_in_MSA_seq2: "+str(pos_in_MSA_seq1)+" , "+str(pos_in_MSA_seq2))
        seq1_from_MSA=[]
        seq2_from_MSA=[]
        j=0
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][pos_in_MSA_seq1]
            seq1_from_MSA.append(found)
            j+=1
        j=0
        while j<=len(MSA_list[0])-1:
            found=MSA_list[0][j][pos_in_MSA_seq2]
            seq2_from_MSA.append(found)
            j+=1
        union=[]
        k=0
        len_max=max(len(seq1_from_MSA),len(seq1_from_MSA))
        while k<=(len_max-1):
            el1=seq1_from_MSA[k]
            el2=seq2_from_MSA[k]
            tuple_like_zip=[el1,el2]
            union.append(tuple_like_zip)
            k+=1
        print("union of the two after merge looks like: "+str(union))
        cost_after_MSA=compute_cost(union,score_matrix,gap_cost)
        if cost_after_MSA==matrix[int(seq1_nr)][int(seq2_nr)]:
            print("integrity test 2 passed for: "+str(seq1_nr)+" and "+ str(seq2_nr))
        else:
            print("Yikes, integrity check 2 did not pas for: "+str(seq1_nr)+" and "+ str(seq2_nr))
            print("Costs were before and after:"+str(matrix[int(node1)][int(node2)])+" and "+str(cost_after_MSA))
            cost_for_suppesed_to_have_been=linear_C(gap_cost,score_matrix,seqs[int(seq1_nr)],seqs[int(seq2_nr)])
            alignment1_str,alignment2_str=linear_backtrack(seqs[int(seq1_nr)], seqs[int(seq2_nr)],cost_for_suppesed_to_have_been, score_matrix, gap_cost)
            alignment1, alignment2 = [*alignment1_str], [*alignment2_str]
            should_have_been= [list(e) for e in zip(alignment1,alignment2)]
            print("should have been:"+ str(should_have_been))
            all_gaps_cols_removed_from_union=[sublist for sublist in union if not all(item == '-' for item in sublist)]
            print(all_gaps_cols_removed_from_union)
            h=0
            while h<=len(should_have_been)-1:
                if should_have_been[h]==all_gaps_cols_removed_from_union[h]:
                    h+=1
                else:
                    print("index of first error: " + str(h) + " out of approximately " + str(len(should_have_been)) + ". The cols are these (should have been, are): " + str(should_have_been[h]) + " and " + str(all_gaps_cols_removed_from_union[h]))
                    sys.exit()
                    h+=1
    total_cost = compute_cost(MSA_list[0], score_matrix, gap_cost)
    print(total_cost)
    return(matrix,min_span_edges_res,in_which_MSA_is_it,MSA_list)



