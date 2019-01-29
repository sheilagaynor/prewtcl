myUpdateCs <-
function (x, K, ws, nstart = 10) 
{
    Cs <- kmeans(x, centers = K, nstart = nstart)$cluster
    sparcl:::UpdateCs(x, K, ws, Cs)
}
