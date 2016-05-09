# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

bp2int = function (path){
  read.str <- paste(readChar(path, file.info(path)$size), collapse=" ")
  nRead <- nchar(read.str)
  readInt <- rep(0, length(read.str))
  for (i in 1:nRead){
    c <- substr(read.str, i, i)
    if (c == 'A'){
      readInt[i] <- 1
    }
    if (c == 'T'){
      readInt[i] <- 2
    }
    if (c == 'C'){
      readInt[i] <- 3
    }
    if (c == 'G'){
      readInt[i] <- 4
    }
  }
  return (list(readInt, nRead))
}