import sys
import gzip
import argparse

correct_char=[False for i in range(ord('z'))]
correct_char[ord('A')]=True
correct_char[ord('C')]=True
correct_char[ord('G')]=True
correct_char[ord('T')]=True
correct_char[ord('U')]=True

def correct_sequence(s):
       for l in s: 
              if not correct_char[ord(l)]: return False
       return True
    
def convert_sequence(s):
       return s.replace("C","T")

# Index all sequences and comments for the bank. For each kmer (read 5'-> 3') store a couple (id of the bank sequence, position in this sequence)
# returns the three data structures (sequences, comments, kmer index)
def index_bank(bankfile, k, convert):
       bank=open(bankfile,"r")
       sequences=[]
       comments=[]
       kmers={}
       id=0
       
       combank=""
       sequence=""
       while True:
              line = bank.readline()
              if not line: break
              if line[0]=='>': # fasta
                     combank=line.rstrip()
                     sequence = bank.readline().rstrip()
              if line[0]=='@': # fastq
                     combank=line.rstrip()
                     sequence = bank.readline().rstrip()
                     bank.readline() # second comment
                     bank.readline() # quality
              if not correct_sequence(sequence): continue
              adapted_sequence = sequence.replace("U","T")
              if convert : adapted_sequence = convert_sequence(adapted_sequence)
             

              
              comments.append(combank)
              sequences.append(sequence)
              
              for i in range (0,len(adapted_sequence)-k):
                  kmer=adapted_sequence[i:i+k]
                  if not kmer in kmers: kmers[kmer]=[]
                  kmers[kmer].append((id,i))

                
              id+=1

       bank.close()
       return sequences,comments,kmers
       


#U (query) versus U bank = not methylated
#C (query) versus U bank (as all 'C' from bank were transformed into U) = methylated

# l1: from query, l2: from bank
# if l1=l2: return 0
# if l1=U and l2=C: return 0:
# else: return 1
def distanceNucleotides (q,b,convert):
       # if q=='T': q='U'# bank = ACGU, query = ACGT
       if q==b: return 0
       if convert and q=='T' and b=='C': return 0 # methylated position (once traitement "bisulfite de sodium")
       return 1
       
# used to print a formated view of a match
def symbolMatchNucleotides (Q,iq,B,ib,convert):

       q=Q[iq]
       b=B[ib]
       if q=='T': q='U'# bank = ACGU, query = ACGT
       if not convert :
           if q==b: return '.'
           else: return ' '
           
           
       if b=='C': b='U'# miRBase 21 convertie in silico (C --> T)
       
       if iq==len(Q)-1 or ib==len(B)-1: 
              if q==b: return '.'
              else: return ' '
       
       

       q_po = Q[iq+1]
       # b_po = B[ib+1]
       if q_po=='T': q_po='U'# bank = ACGU, query = ACGT
       # if b_po=='C': b_po='U'# miRBase 21 convertie in silico (C --> T)
       
       if q=='C':
              if q_po=='G': return 'Y'
              else: return 'C'
              
       
       
       if q==b: return '.'
       

       return ' '
       
       



# compare two (same length) sequences Q and B. At least "minimal_range_query" of the query must match
def compareQueryAndRef (Q,B,start_on_bank,convert,minimal_range_query,threshold=0):
       dist=0
       start=0
       minimal_span = minimal_range_query*len(Q)//100
       maximal_start = (100-minimal_range_query)*(len(Q)+threshold)//100
       minimal_start = -maximal_start
       
       if start_on_bank<minimal_start: return threshold+1 # cannot align a sufficient long part of the query (starting before B)
       if (len(B)-start_on_bank)<minimal_span: return threshold+1 # cannot align a sufficient long part of the query (ending after B)
       if start_on_bank<0:
           
           B_tmp=["X" for i in range(-start_on_bank)]
           B_tmp+=B[0:len(Q)+start_on_bank]
           B=B_tmp
       else: 
           B=B[start_on_bank:start_on_bank+len(Q)]
       
       #Starting from begining 
       for i in range(min(len(Q), len(B))):
              dist+=distanceNucleotides(Q[i],B[i],convert)
              if dist>threshold: 
                     if (i-start)>minimal_span: 
                            return dist-1
                     if i>maximal_start: return dist # no way we can obtain at least 90% aligned
                     
                     # else: start back from i
                     start=i+1 
                     dist=0
       if min(len(Q), len(B))-start<minimal_span: return dist+1
       return dist
       
def printAMatch(B,comquery,Q,start, convert):
    #ADDED 12/12/2016
       if start<0: 
           return #Finaly not interested by matching starting before the premi
       if start+len(Q)>len(B):
           return #Finaly not interested by matching starting after the premi
    #END ADDITION 12/12/2016
       spaces=""
       if start >=0:
           for i in range(start): spaces+=" "
       else:
           Q=Q[-start:]
           start=0
       print (spaces+Q,"\t",comquery)
       print (spaces, end="")
       for i in range(min(len(Q),len(B)-start)): 
              print (symbolMatchNucleotides(Q,i,B,i+start,convert),end="")
       print ()
       


       

# Compare a full query sequence to the whole indexed bank. Returns the best alignement (id of the reference and position on the reference)
def query(sequences,comments,kmers, query,k,convert,minimal_range_query,threshold=0):
       # print (query)
       query = query.replace("U","T")
       if convert : adapted_query = convert_sequence(query) # Used only for finding seeds. 
       else: adapted_query = query
       best_result_ids=[]# no reference at no position
       best_result_value = threshold+1 # with a bad score (threshold+100)
       ref_tested={} # stores the key (bank id) value (positions) of the already tested mapping positions
       for kmer_on_query in range (0,len(query)-k+1):
              kmer=adapted_query[kmer_on_query:kmer_on_query+k] # kmers (seeds) are taken from converted sequences (if convert==True)
              if not kmer in kmers: continue
              for (bank_sequence_id,kmer_on_bank) in kmers[kmer]:
                     # ------------------- bank
                     #    ---------------- query
                     #        <kmer>: (kmer_on_query= 4, kmer_on_bank=7 --> start_on_bank=3)
                     start_on_bank=kmer_on_bank-kmer_on_query
                     
                     # if start_on_bank<0 : continue
                     
                     
                     if bank_sequence_id in ref_tested and start_on_bank in ref_tested[bank_sequence_id]: continue # Already tested
                     
                     if not bank_sequence_id in ref_tested: ref_tested[bank_sequence_id]=[]
                     ref_tested[bank_sequence_id].append(start_on_bank) # We mark this couple bank_sequence_id, start_on_bank as tested
                     bank_sequence = sequences[bank_sequence_id].replace("U","T")
                     dist = compareQueryAndRef (query,bank_sequence, start_on_bank,convert,minimal_range_query,threshold) # the distances are made on the non converted sequences (only ACGT-based and non ACGU-based)
                     # print ("dist = "+str(dist))
                     if dist>threshold: continue

                     if (dist==best_result_value): 
                            best_result_ids.append((bank_sequence_id,start_on_bank))
                     if (dist<best_result_value): 
                            best_result_value=dist
                            best_result_ids=[]
                            best_result_ids.append((bank_sequence_id,start_on_bank))
                            # print ("hey", start_on_bank)
 #                            for i in range(kmer_on_query): print(" ",end="")
 #                            print (kmer)
 #                            print (query)
 #                            print (bank_sequence[start_on_bank:start_on_bank+len(query)])
                            # print (bank_sequence)
       return best_result_ids
              
       
              
def compare_all_queries(queryfile,sequences,comments,kmers,k,convert,minimal_range_query,threshold=0):
       nb_queries=0
       matches={}
       
       queries=open(queryfile,"r")
       comquery=""
       sequence=""
       while True:
              nb_queries+=1
              if(nb_queries%10000==0): print (nb_queries, "queries treated")
              line = str(queries.readline()).lstrip('b\'').rstrip()
              if not line: break
              if line[0]=='>': # fasta
                     comquery=line
                     sequence = queries.readline().rstrip()
              if line[0]=='@': # fastq
                     comquery=line.rstrip()
                     sequence = str(queries.readline()).rstrip().lstrip('b\'')
                     queries.readline() # second comment
                     queries.readline() # quality

              if correct_sequence(sequence):
                     couples_bank_sequence_id_bank_sequence_position=query(sequences,comments,kmers, sequence,k,convert,minimal_range_query,threshold)
                     if len(couples_bank_sequence_id_bank_sequence_position)==0: continue
                     for (bank_sequence_id,bank_sequence_position) in couples_bank_sequence_id_bank_sequence_position:
                            if bank_sequence_id not in matches: matches[bank_sequence_id]=[]
                            matches[bank_sequence_id].append((comquery,sequence,bank_sequence_position))                     #
                            # seqbank = sequences[bank_sequence_id]
                            # combank = comments[bank_sequence_id]
                            # print (seqbank+"\t"+combank+"\t",bank_sequence_position)
                            # printAMatch(seqbank,comquery,sequence,bank_sequence_position)
                     
       queries.close()
       return matches

def print_results(sequences,comments,matches, convert):
       for bank_sequence_id in matches: 
                     seqbank = sequences[bank_sequence_id]
                     combank = comments[bank_sequence_id]
                     print (seqbank+"\t"+combank)
                     for (comquery,sequence,bank_sequence_position) in matches[bank_sequence_id]:
                            printAMatch(seqbank,comquery,sequence,bank_sequence_position, convert)
                      
def main():
    parser = argparse.ArgumentParser(description="Maps short sequences to a reference bank. If required,  \'T\'s from queries match \'C\'s from the bank.  ")
    parser.add_argument("input_bank_file", type=str,
                        help="input fasta or fastq bank file" )
    parser.add_argument("input_query_file", type=str,
                        help="input fasta or fastq query file" )
    parser.add_argument("converted", type=str,
                        help="chose \"True\" or \"False\". False: usual mapping. True:  \"T\"s from queries match \"C\"s from the bank.")
    parser.add_argument("-k", type=int, dest='k',
                        help="kmer size [default: 12]", default=12 )
    parser.add_argument("-t", type=int, dest='t',
                        help="Maximal number authorized substitution [default: 0]", default=0 )
    parser.add_argument("-span", type=int, dest='s',
                        help="The portion of a read mapped on a reference may be lower than 100 percent. Span (in 0-100) provides this minimal percentage value [Default 90]", default=90)
    
    args = parser.parse_args()
    k=args.k
    t=args.t
    s=args.s
    if args.converted=="True" or args.converted=="true" or args.converted=="T" or args.converted=="t":
        convert=True
    else:
        convert=False
    
    print("indexation")
    sequences,comments,kmers=index_bank(args.input_bank_file, k, convert)
    print("querying, with threshold="+str(t))
    matches=compare_all_queries(args.input_query_file,sequences,comments,kmers,k,convert,s,t)
    print ("\n\n\t\t ******** RESULTS ********\n\n")
    print_results(sequences,comments,matches,convert)
    
if __name__ == "__main__":
    main()    
