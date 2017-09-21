'''
Joaquim Claivaz
170526

Inference if each gene of ortholog pair have at least one domain inferred by pfamscan
'''


#Ubuntu
setwd('/home/jclaivaz/Bureau')
#Windows
setwd('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference/other/')

#DROME_DROPS
data_tmp_DROME_DROPS=read.table('status_output_DROME_DROPS')

for (i in 1:dim(data_tmp_DROME_DROPS)[1]){
  if (data_tmp_DROME_DROPS$V2[i]=='yes'){
    if (data_tmp_DROME_DROPS$V4[i]=='yes'){
      data_tmp_DROME_DROPS$status[i]='both'
    }else{
      data_tmp_DROME_DROPS$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROPS$V4[i]=='yes'){
      data_tmp_DROME_DROPS$status[i]='partial'
    }else{
      data_tmp_DROME_DROPS$status[i]='none'
    }
  }
}

write.csv(data_tmp_DROME_DROPS, 'DROME_DROPS_status')

table(data_tmp_DROME_DROPS$status)

#DROME_DROYA
data_tmp_DROME_DROYA=read.table('status_output_DROME_DROYA')

for (i in 1:dim(data_tmp_DROME_DROYA)[1]){
  if (data_tmp_DROME_DROYA$V2[i]=='yes'){
    if (data_tmp_DROME_DROYA$V4[i]=='yes'){
      data_tmp_DROME_DROYA$status[i]='both'
    }else{
      data_tmp_DROME_DROYA$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROYA$V4[i]=='yes'){
      data_tmp_DROME_DROYA$status[i]='partial'
    }else{
      data_tmp_DROME_DROYA$status[i]='none'
    }
  }
}
write.csv(data_tmp_DROME_DROYA, 'DROME_DROYA_status')

table(data_tmp_DROME_DROYA$status)

#DROME_DROAN
data_tmp_DROME_DROAN=read.table('status_output_DROME_DROAN')

for (i in 1:dim(data_tmp_DROME_DROAN)[1]){
  if (data_tmp_DROME_DROAN$V2[i]=='yes'){
    if (data_tmp_DROME_DROAN$V4[i]=='yes'){
      data_tmp_DROME_DROAN$status[i]='both'
    }else{
      data_tmp_DROME_DROAN$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROAN$V4[i]=='yes'){
      data_tmp_DROME_DROAN$status[i]='partial'
    }else{
      data_tmp_DROME_DROAN$status[i]='none'
    }
  }
}

write.csv(data_tmp_DROME_DROAN, 'DROME_DROAN_status')

table(data_tmp_DROME_DROAN$status)

#DROME_DROMO
data_tmp_DROME_DROMO=read.table('status_output_DROME_DROMO')

for (i in 1:dim(data_tmp_DROME_DROMO)[1]){
  if (data_tmp_DROME_DROMO$V2[i]=='yes'){
    if (data_tmp_DROME_DROMO$V4[i]=='yes'){
      data_tmp_DROME_DROMO$status[i]='both'
    }else{
      data_tmp_DROME_DROMO$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROMO$V4[i]=='yes'){
      data_tmp_DROME_DROMO$status[i]='partial'
    }else{
      data_tmp_DROME_DROMO$status[i]='none'
    }
  }
}

write.csv(data_tmp_DROME_DROMO, 'DROME_DROMO_status')

table(data_tmp_DROME_DROMO$status)

#DROME_DROSI
data_tmp_DROME_DROSI=read.table('status_output_DROME_DROSI')

for (i in 1:dim(data_tmp_DROME_DROSI)[1]){
  if (data_tmp_DROME_DROSI$V2[i]=='yes'){
    if (data_tmp_DROME_DROSI$V4[i]=='yes'){
      data_tmp_DROME_DROSI$status[i]='both'
    }else{
      data_tmp_DROME_DROSI$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROSI$V4[i]=='yes'){
      data_tmp_DROME_DROSI$status[i]='partial'
    }else{
      data_tmp_DROME_DROSI$status[i]='none'
    }
  }
}

write.csv(data_tmp_DROME_DROSI, 'DROME_DROSI_status')

table(data_tmp_DROME_DROSI$status)

#DROME_DROVI
data_tmp_DROME_DROVI=read.table('status_output_DROME_DROVI')

for (i in 1:dim(data_tmp_DROME_DROVI)[1]){
  if (data_tmp_DROME_DROVI$V2[i]=='yes'){
    if (data_tmp_DROME_DROVI$V4[i]=='yes'){
      data_tmp_DROME_DROVI$status[i]='both'
    }else{
      data_tmp_DROME_DROVI$status[i]='partial'
    }
  }else{
    if (data_tmp_DROME_DROVI$V4[i]=='yes'){
      data_tmp_DROME_DROVI$status[i]='partial'
    }else{
      data_tmp_DROME_DROVI$status[i]='none'
    }
  }
}

write.csv(data_tmp_DROME_DROVI, 'DROME_DROVI_status')

table(data_tmp_DROME_DROVI$status)
