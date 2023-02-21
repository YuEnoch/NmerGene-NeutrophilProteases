#YuEnoch 20-02-2023
#MakeNNKCombinations.py

#Purpose: creates all possible n-mer combinations
#         in this case, a 6-mer (but add for loops manually)

AAList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y']

output = open("NNK6_combinations.txt", "w")
for i in range(len(AAList)):
    for j in range(len(AAList)):
        for k in range(len(AAList)):        
            for l in range(len(AAList)):    
                for m in range(len(AAList)):    
                    for n in range(len(AAList)):    
                        str = AAList[i] + AAList[j] + AAList[k] + AAList[l] + AAList[m] + AAList[n]
                        print(str, file = output)
output.close()
