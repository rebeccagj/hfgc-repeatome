# published here, by these authors, got file with this wget command:
# wget [remember to add command here]

setwd("/Users/Rebecca/Downloads/") #set wd to directory with dros txt file outputs

fxnRF = function(obj_arg, file_arg){
  obj_arg = read.delim(file = file_arg, header = T) #tab delim
}

# whole ####
whole1 = fxnRF(whole1, "GSM2650779_whole1Genes.txt")
whole2 = fxnRF(whole2, "GSM2650780_whole2Genes.txt")
whole3 = fxnRF(whole3, "GSM2650781_whole3Genes.txt")
whole4 = fxnRF(whole4, "GSM2650782_whole4Genes.txt")
whole5 = fxnRF(whole5, "GSM2650783_whole5Genes.txt")
#can do a merge later as independant practice


# c4da ####
c4da1 = fxnRF(c4da1, "GSM2650775_c4da1Genes.txt")
c4da2 = fxnRF(c4da2, "GSM2650776_c4da2Genes.txt")
c4da3 = fxnRF(c4da3, "GSM2650777_c4da3Genes.txt")
c4da4 = fxnRF(c4da4, "GSM2650778_c4da4Genes.txt")

allc4da = merge.data.frame(c4da1, c4da2, all.x = T, all.y = T)
allc4da = merge.data.frame(allc4da, c4da3, all.x = T, all.y = T)
allc4da = merge.data.frame(allc4da, c4da4, all.x = T, all.y = T)

write.csv(allc4da, file = "20180308_allc4da.csv") #write merged data.frame to csv file
allc4da = read.csv("20180308_allc4da.csv", row.names = 1) #row names var to prevent R from appending X column of numbers