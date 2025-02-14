#Bioinformatics pipeline:
#I used this for human proteins but should work with other databases and stuff assuming you know what to change
#run this line to run the program:
#protname [query] -blastcom [blastp/n] -musclecom muscle.exe -fastafile [database containing query] -databasefasta [blast search database]

from Bio import SeqIO
import random
import os
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import *
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
#Function to return command line arguments and print help message
def get_arguments_a3():
    import argparse
    parser = argparse.ArgumentParser(description="Assignment 3")
    parser.add_argument("-protname", type=str, help="Protein name",required=True)
    parser.add_argument("-blastcom", type=str, help="Blast executable",required=True)
    parser.add_argument("-musclecom", type=str, help="muscle executable",required=True)
    parser.add_argument("-fastafile", type=str, help="Source fasta file",required=True)
    parser.add_argument("-databasefasta", type=str, help="Target database file",required=True)
    args = parser.parse_args()
    return [arg for arg in vars(args).values()]


print("To add mutiple proteins, sepearate each protein with a comma")
if __name__ == "__main__":
    input_args = get_arguments_a3()
    protein_name = input_args[0].split(",") #the code here is changed so as to input multiple proteins and seperating them based on a comma(AI USED)
    blast_exe = input_args[1]
    muscle_exe = input_args[2]
    human_fasta = input_args[3]
    database_fasta = input_args[4]
    
    print("Pipeline running with protein name",protein_name)
    print("Running with BLAST as",blast_exe)
    print("Running with muscle as",muscle_exe)
    print("Source fasta:",human_fasta)
    print("Database to search against:",database_fasta)








    for pn in protein_name:    #The for loop goes through all the proteins specified in the input and performs the entire analysis for each protein(AI USED)
        rlist=[]    #List for storing the records from database_fasta where all the human proteins are stored
        with open(human_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                rlist.append(record.id)    #List of all human protein sequences from the database


        RND=random.randint(0,10000)    #Generates a random number to use as file names to differentiate between the outputs

        ###BLAST PART###

        blast_input=f"temp_{RND}_blast_input.txt"    #input file for the blast search
        blast_output=f"temp_{RND}_blast_output.tsv"  #output file for the blast search

        with open(blast_input,"w") as f:    #opens the blast input file and "w" indicates that the file will be written
            for record in SeqIO.parse(human_fasta,"fasta"):     #goes through all sequences in the file
                name=record.id.split("|")[2]    #splits the name of the human protein to make it cleaner in the output, splits after the "|" are accessed by [2] (e.g:>sp|Q86UQ4|ABCAD_HUMAN becomes ABCAD_HUMAN)  )(AI Used)
                if name==pn:    #if the protein in the database matches the one 
                    f.write(f">{record.id}\n{record.seq}\n")    #the Id and the sequence is written the blast input file

        command = blast_exe + " -db UP_mammals -query "+ blast_input +" -out " + blast_output + " -outfmt 6"    #command that runs the BLAST Search 
        os.system(command)
        f.close()
        
        #(AI used for lines 75-79, modified)
        blastdf=pd.read_csv(blast_output,sep="\t", header=None)    #using Pandas to sort the output file into a tab-seperated file
        s_col=blastdf.iloc[:,[0,1,10,11]]    #Isolates the column for Input Sequence, Hit Sequence, E Value and BLAST Scores
        s_col.columns=["Input Sequence", "Hit_Sequence","E_value","BLAST Score"]   #Gives names to the correspoding columns
        filtered_s_col=s_col[s_col["E_value"] <= 1e-5]# Gives the entire row that has an E-Value of less than 1E-5
        filtered_s_col.to_csv(f'filtered_blast_output_{RND}.tsv', sep='\t', index=False)    #Converts the file to CSV 






        id_to_record = {record.id: record for record in SeqIO.parse(database_fasta, "fasta")}    #this dictonary stores all the record ids and sequences from the database fasta as keys and values respectively(AI USED)

        species_name=[]    #list to store the names after removing the "_"
        for record_id in id_to_record.keys():     #goes through all the keys in the dictonary
            species_name.append(record_id.split("_")[-1])    #each species is seprated based on "_" and then the next part is taken using [-1] which is the name itself (eg: >tr|A0A3B2WB67|A0A3B2WB67_MOUSE to MOUSE)
        
        records_hit=[]    #list that stores all the id of the sequences that are also present in the blast output
        seen_species=set()    #this is a special type of list that only adds elements that do not repeat, we do this because there are multiple ids with the same name in the blast search, it only stores the first element it encounters and discards its repeats(AI USED)

        for _, row in filtered_s_col.iterrows():    #goes through all the rows in the blast search(AI USED)
            hit_id=row["Hit_Sequence"]   #takes the value in the column hit sequence for that row
            species=hit_id.split("_")[-1]#splits the value obtained based on "_" to get only species name(eg: >tr|A0A3B2WB67|A0A3B2WB67_MOUSE to MOUSE)
            if species not in seen_species and hit_id in id_to_record: #checks if the blast search name is not present in our  seen_species list and also checks if the hit id from the blast search is in the database fasta(AI USED)
                records_hit.append(id_to_record[hit_id])
                seen_species.add(species)   #add the species to our unique list 






        ###Multiple Sequence Alignment Part###
        msa_input=f"temp_{RND}_msa_input.fasta" #file stores the result from our blast part
        SeqIO.write(records_hit, msa_input,"fasta")   
        msa_output=f"temp_{RND}_msa_output.fasta" #files stores results from the msa done by muscle_exe
        command_msa = muscle_exe + " -align " + msa_input + " -output " + msa_output #command that runs alignment done by muscle
        os.system(command_msa)

        #AI USED from line 116-139, modified
        from collections import Counter  #this library contains counter that counts the number of elements, it is unique to a regular counter cause it gives a dictonary with the element and its frequency (eg:{'a': 3, 'b': 2, 'c': 1})
        al=AlignIO.read(msa_output, "fasta")# reads the aligned sequneces obtained from the muscle output
        #con=[]#contains the frequency of amino acids appeard in the same position across multiple sequences
        con_scores=[]# calculates the conservation score i.e the frequency of the amino acid appearing in the same position across multiple sequences
        aminoacids=[]
        con_dict={}# dictonary used to store the amino acid as a key and their total frequency across all positions, this is for the pie chart
        for i in range(0, len(al[0])):#For loop that goes through all the amino acids in the seqeunces(all are of the same lenght, after alignment)
            col = [record[i] for record in al]#takes the amino acid from all the seqeunces at the same position, (col is a list of those amino acids), it also contains gaps(-)
            col=[aa for aa in col if aa!="-"]#goes through all amino acids in the col list and proceeds if the amino acid is not a gap
            count = Counter(col)#this dictonary uses counter to which has the amino acid as key and its freqency as value of that key
            mcr = count.most_common(1)[0][1]#nost common is a function that returns element with the most freqeuncy, here (1), indicates the most freqeuent element and its frequency, [0] converts it to a list that contains both the amino acid and the freqeuncy and [1] gets only the frequency
            mcaa=count.most_common(1)[0][0]#same principle as the above but this takes the amino acid instead
            aminoacids.append(mcaa)#list that stores all the amino acid
            con_score=mcr/len(al)#calculates conservation score of the amino acid by dividing the frequency by the enitre length of the sequence
            con_scores.append(con_score)#list that stores the scores calculated from above
            con_dict[mcaa]=con_dict.get(mcaa,0)+con_score#gets the conservation score from the con_dict dictonary, if there element is present in the dictornary add the con_score, if not add 0, this calculates the total conservation score of the amino acid across all sequences, this is to get the most common amino acid for the pie chart
        sorted_con_dict=sorted(con_dict.items(),key=lambda x: x[1],reverse=True)#sorts all the items in the dictonary in descending order (reverse=True), key=lambda x:x[1] gets only the value i.e conservation score for the amino acid, items gets the key
        top5=sorted_con_dict[:5]#takes the first five elements of the sorted_con_dict, this includes the key and the value
        labels = [item[0] for item in top5]#list that contains the amino acids, to act as labels in the pie plot
        sizes = [item[1] for item in top5]#list that contains the correspoding conservation scores

        con_plot_pie=plt.figure(figsize=(6, 6))#plots a pie chart containing the top 5 most common amino acids across all sequences
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=plt.cm.Paired.colors)#sizes is the conservation score, labels are the correspoding amino acid, autopct shows the percentage of the slice in the pie chart, colors gives color
        plt.title(f"Top 5 Most Conserved Amino Acids in Alignment", fontsize=14, fontweight="bold")




    
        ###Consensus Plot###
        x=range(1, len(con_scores)+1)#x axis that contains the all the positions of the amino acids and 1 is added cause it starts from 1 (0 doesnt make sense to start with when referring to amino acid postion)
        y=con_scores#conservation score of each amino acid
        max_width=50# the maximum width the plot can take in inches
        fig_width = min(len(con_scores) * 0.1, max_width)#(AI USED)
        con_scores_plot=plt.figure(figsize=(fig_width, 6))#the actual line plot 
        plt.plot(x, y, color="blue", linewidth=2)#parameters of the line 
        plt.title(f"Conservation of Amino Acids",  fontsize=14, fontweight="bold")
        plt.xlabel("Amino Acid Position")#label for x axis
        plt.ylabel("Conservation")#label for y axis
        plt.title(f"Consensus Plot for Alignment", fontsize=14, fontweight="bold")#plot title
        plt.grid(False)#removes grid from the plot(AI USED)
        plt.tight_layout()#ensures that all elements fit within the figure size and prevents overlap of elements(AI USED)







        ###Phylogenetics Part##
        calculator=DistanceCalculator("blosum62")#used to get the distance matrix calculated using blosum62 method
        dm=calculator.get_distance(al)#gets the distance value from the matrix
        constructor=DistanceTreeConstructor()#method for calulatiing the tree
        tree=constructor.nj(dm)#contains the tree object calculated using the neighbour joining method
        matplotlib.rc('font', size=5)#rc allows for modification of the standard matplot, this is used to construct the tree
        fig,ax = plt.subplots(figsize=(6,6))#makes the plot
        nj_plot=fig#called later

        for clade in tree.find_clades():#for loop that is used to access the clade terminals and the branches
            if clade.name and "Inner" in clade.name:#takes the inner that is present and replaces it with "" to makes the tree look cleaner
                clade.name = ""
            if clade.is_terminal():#this takes the id that is present at the end of a branch
                mod_name=clade.name.split("|")[-1]#splits the id based on "|" and takes the element that is after that(AI USED)
                mod_name1=mod_name.replace("|","_")#in that element | is replaced by "_"(AI USED)
                clade.name=mod_name1#replaces the clade name with the modified name(AI USED)


            bl=clade.branch_length#gets the lenght of the branch and changes color based on the length(AI USED)
            if bl<0.01:
                color="brown"
            elif bl<0.02:
                color="orange"
            elif bl<0.03:
                color="maroon"
            else:
                color="red"
            clade.color=color#converts the current branch color with the color that is obtained from the above statements
            
        Phylo.draw(tree,do_show=False, axes=ax)#renders the tree, do show, shows the plot immediately after this is executed, this is set to false, it will show at the end, axes=shows the axis
        ax.set_title(f"Phylogenetics Tree--Neighbour Joining", fontsize=14, fontweight="bold")#plot title

        ###Works the same as above, only difference is that the tree is based on UPGMA method rather that Neighbour Joining method for construction
        treeupgma=constructor.upgma(dm)

        matplotlib.rc('font', size=5)
        fig1,ax1 = plt.subplots(figsize=(6,6))
        upgma_plot=fig1

        for clade in treeupgma.find_clades():
            if clade.name and "Inner" in clade.name:
                clade.name = ""
            if clade.is_terminal():
                modu_name=clade.name.split("|")[-1]#AI USED
                modu_name1=modu_name.replace("|","_")#AI USED
                clade.name=modu_name1

            
            bl=clade.branch_length#AI USED
            if bl<0.01:
                color="brown"
            elif bl<0.02:
                color="orange"
            elif bl<0.03:
                color="maroon"
            else:
                color="red"
            clade.color=color

        Phylo.draw(treeupgma,do_show=False, axes=ax1)
        ax1.set_title(f"Phylogenetics Tree--UPGMA",  fontsize=14, fontweight="bold")






        #Additional Things###


        ##Histogram of the BLAST scores##
        #shows how many of the results have that blast score
        BScores=filtered_s_col["BLAST Score"]#takes the column containing all the blast scores
        BScores=BScores.tolist()#converts the column into a list
        bsp=plt.figure(figsize=(6,6))#makes the figure
        plt.hist(BScores, bins=20, color="black")#adds the histogram(AI USED)
        plt.xlabel("BLAST Score")
        plt.ylabel("Frequency")
        plt.title(f"Histogram of the BLAST Scores",  fontsize=14, fontweight="bold")
        plt.grid(True)#shows the grid
        plt.xticks(range(int(min(BScores)), int(max(BScores)) + 1, int((max(BScores)-min(BScores))//10)))#calculates the distances between the points on the x axis, takes teh minimum blast score and the highest blast score and calculates the appropriate range between them(AI USED)

    

        from matplotlib.backends.backend_pdf import PdfPages#library used to make the multipage pdf(AI USED)

        #Table of top 5 blast scores

        sorted_blast_scores=filtered_s_col.sort_values(by="BLAST Score", ascending=False)#takes the blast search results are arranges them accoring to blast score, from highest to lowest(AI USED)
        top_5_blast_scores=sorted_blast_scores.head(5)#takes the first five items in the table
        top_5_blast_scores_df=top_5_blast_scores[["Hit_Sequence", "E_value", "BLAST Score"]]#gives titles to the 5 items
        num_cols = len(top_5_blast_scores_df.columns)
        fig_width = max(5, num_cols * 3)#(AI USED)
        fig, ax=plt.subplots(figsize=(fig_width,3))#figure for the table
        ax.axis('off')#removes the axes
        ax.set_title(f"Top 5 Blast Hits for {pn}",fontsize=14, fontweight="bold")#Title of the table
        table = ax.table(cellText=top_5_blast_scores_df.values, colLabels=top_5_blast_scores_df.columns, loc='center')#the table with columns being the id, celltext being the blast scores, and loc is the position of the table on the page
        table.auto_set_font_size(False)#all the text remains the same size in the table(AI USED)
        table.set_fontsize(10)#(AI USED)
        table.auto_set_column_width(col=list(range(len(top_5_blast_scores_df.columns))))#sets the column width according to text size, col is the columnns where the size is adjusted, in this case it is all the entries starting from the 0th row
        blast_table_plot=fig#contains the figure


        #Intro Page that contains the title, input protein, name of the database and programs used
        fig,ax=plt.subplots(figsize=(6,6))#AI used
        ax.axis("off")
        title=f"{pn} Analysis Report"#title
        ax.text(0.5, 1, title, ha="center", fontsize=20, fontweight="bold")#centers the title based on the coordinates(x, y), ha= horizontal alignment
        ax.text(0.2,0.8, "Input Protein", ha="center", fontsize=12)#on the left side
        ax.text(0.8,0.8, f"{pn}", ha="center", fontsize=12)#protein name on the right
        ax.text(0.2,0.7, "Searching in", ha="center", fontsize=12)# on the left
        ax.text(0.8,0.7,f"{database_fasta}", ha="center", fontsize=12)#database name on the right
        ax.text(0.2,0.6,"Programs used", ha="center", fontsize=12)#on the left
        ax.text(0.8,0.6, f"{blast_exe}, {muscle_exe}", ha="center", fontsize=12)#programs used
        intro=fig




        fig,ax=plt.subplots(figsize=(6,6))
        ax.axis("off")
        title=f"Raw Files"#title
        ax.text(0.5, 1, title, ha="center", fontsize=20, fontweight="bold")#centers the title based on the coordinates(x, y), ha= horizontal alignment
        ax.text(0.2,0.8, "Blast Input", ha="center", fontsize=8)#on the left side
        ax.text(0.8,0.8, f"temp_{RND}_blast_input.txt", ha="center", fontsize=8)#blast input name on the right
        ax.text(0.2,0.7, "Blast Output", ha="center", fontsize=8)#on the left side
        ax.text(0.8,0.7, f"temp_{RND}_blast_output.txt", ha="center", fontsize=8)#blast output name on the right
        ax.text(0.2,0.6, "Filtered Blast Output(contains E-values of less that 1E-5)", ha="center", fontsize=8)#on the left side
        ax.text(0.8,0.6, f"filtered_blast_output_{RND}.tsv", ha="center", fontsize=8)#filterd blast output
        ax.text(0.2,0.5, "Alignment Input", ha="center", fontsize=8)#on the left side
        ax.text(0.8,0.5, f"temp_{RND}_msa_input.fasta", ha="center", fontsize=8)#msa input name on the right
        ax.text(0.2,0.4, "Alignment Output", ha="center", fontsize=8)#on the left side
        ax.text(0.8,0.4, f"temp_{RND}_msa_output.fasta", ha="center", fontsize=8)#msa output name on the right
        lastpage=fig

        #compiling into a multipage pdf
        with PdfPages(f"Analysis Report for {pn}.pdf") as pdf:#name of the final pdf
            pdf.savefig(intro)#intro page
            plt.close(intro)#closed to prevent memory usage

            pdf.savefig(blast_table_plot)#table of top 5 blast results
            plt.close(blast_table_plot)

            pdf.savefig(bsp)#histogram of blast results
            plt.close(bsp)

            pdf.savefig(con_scores_plot)#consensus plot
            plt.close(con_scores_plot)

            pdf.savefig(con_plot_pie)#pie plot of top 5 amino acids
            plt.close(con_plot_pie)

            pdf.savefig(nj_plot)#nj tree
            plt.close(nj_plot)

            pdf.savefig(upgma_plot)#upgma tree
            plt.close(upgma_plot)

            pdf.savefig(lastpage)#last page
            plt.close(lastpage)

       





    







    


        

    




    
        
        




