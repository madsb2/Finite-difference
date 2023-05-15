#2.9a
alpha <- function(i,pos,dx){
  if(pos==-2){
    return (dx[i] / (dx[i - 1] * (dx[i - 1] + dx[i])))
  }else if(pos==-1){
    return ((-dx[i - 1] - dx[i]) / (dx[i - 1] * dx[i]))
  }else if(pos==0){
    return ((dx[i - 1] + 2 * dx[i]) / (dx[i] * (dx[i - 1] + dx[i])))
  }else{
    return(NA)
  }
}

#2.9b
beta <- function(i, pos, dx){
  if(pos==-1){
    return (-dx[i + 1] / (dx[i] * (dx[i] + dx[i + 1])))
  }else if(pos==0){
    return ((dx[i + 1] - dx[i]) / (dx[i] * dx[i + 1]))
  }else if(pos==1){
    return (dx[i] / (dx[i + 1] * (dx[i] + dx[i + 1])))
  }else{
    return(NA)
  }
}

#2.9c
gamma <- function(i, pos, dx){
  if(pos==0){
    return ((-2 * dx[i + 1] - dx[i + 2]) / (dx[i + 1] * (dx[i + 1] + dx[i + 2])))
  }else if(pos==1){
    return ((dx[i + 1] + dx[i + 2]) / (dx[i + 1] * dx[i + 2]))
  }else if(pos==2){
    return (-dx[i + 1] / (dx[i + 2] * (dx[i + 1] + dx[i + 2])))
  }else{
    return(NA)
  }
}

#2.10
delta <- function(i, pos, dx){
  if(pos==-1){
    return (2 / (dx[i] * (dx[i] + dx[i + 1])))
  }else if(pos==0){
    return (-2 / (dx[i] * dx[i + 1]))
  }else if(pos==1){
    return (2 / (dx[i + 1] * (dx[i] + dx[i + 1])))
  }else{
    return(NA)
  }
}
