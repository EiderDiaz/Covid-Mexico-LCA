clusterStability <- function(data=NULL, clustermethod=NULL, randomTests = 20, trainFraction = 0.5,pac.thr=0.1, ...)
{
  clusterLabels <- list();
  randomSamples <- list();
  numberofClusters <- 0;
  testCounts <- numeric(nrow(data))
  randomSeeds <- sample(randomTests);

  for (i in 1:randomTests)
  {
    randomSamples[[i]] <- sample(nrow(data),trainFraction*nrow(data));
    mod1 <- clustermethod(data[randomSamples[[i]],],...);
    clusterLabels[[i]] <- predict(mod1,data);
    names(clusterLabels[[i]]$classification) <- rownames(data)
    plot(data[,1:2],col = clusterLabels[[i]]$classification,main=sprintf("%d",i));
    numberofClusters <- numberofClusters + length(table(clusterLabels[[i]]$classification))
    testCounts[-randomSamples[[i]]] <- testCounts[-randomSamples[[i]]] + 1;
    set.seed(randomSeeds[i]);
  }
  numberofClusters <- numberofClusters/randomTests;
  cat("Done Testing:")
  randIndex <- numeric();
  jaccIndex <- numeric();
  meanJaccard <- numeric();
  jaccardpoint <- numeric(nrow(data));
  jaccardpointcount <- numeric(nrow(data));
  trainrandIndex <- numeric();
  trainjaccIndex <- numeric();
  trainmeanJaccard <- numeric();
  trainjaccardpoint <- numeric(nrow(data));
  trainjaccardpointcount <- numeric(nrow(data));
  for (i in 1:(randomTests - 1))
  {
    for (j in (i + 1):randomTests)
    {
      outsamples <- unique(c(randomSamples[[i]],randomSamples[[j]]))
      if ((nrow(data) - length(outsamples)) > 10)
      {
        randIndex <- c(randIndex,adjustedRandIndex(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]));
        jaccard <- jaccardMatrix(clusterLabels[[i]]$classification[-outsamples],clusterLabels[[j]]$classification[-outsamples]);
        jaccIndex <- c(jaccIndex,jaccard$balancedMeanJaccard);
        meanJaccard <- c(meanJaccard,mean(jaccard$elementJaccard));
        jaccardpoint[-outsamples] <- jaccardpoint[-outsamples] + jaccard$elementJaccard;
        jaccardpointcount[-outsamples] <- jaccardpointcount[-outsamples] + 1;
      }
      insamples <- randomSamples[[i]][randomSamples[[i]] %in% randomSamples[[j]]]
      if ((nrow(data) - length(insamples)) > 10)
      {
        trainrandIndex <- c(trainrandIndex,adjustedRandIndex(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]));
        trainjaccard <- jaccardMatrix(clusterLabels[[i]]$classification[insamples],clusterLabels[[j]]$classification[insamples]);
        trainjaccIndex <- c(trainjaccIndex,trainjaccard$balancedMeanJaccard);
        trainmeanJaccard <- c(trainmeanJaccard,mean(trainjaccard$elementJaccard));
        trainjaccardpoint[insamples] <- trainjaccardpoint[insamples] + trainjaccard$elementJaccard;
        trainjaccardpointcount[insamples] <- trainjaccardpointcount[insamples] + 1;
      }
    }
  }
  cat("After Jacckard:")
  jaccardpoint[jaccardpointcount > 0] <- jaccardpoint[jaccardpointcount > 0]/jaccardpointcount[jaccardpointcount > 0];
  names(jaccardpoint) <- rownames(data);
  trainjaccardpoint[trainjaccardpointcount > 0] <- trainjaccardpoint[trainjaccardpointcount > 0]/trainjaccardpointcount[trainjaccardpointcount > 0];
  names(trainjaccardpoint) <- rownames(data);

  testConsesus <- matrix(0,nrow = nrow(data), ncol = nrow(data))
  colnames(testConsesus) <- rownames(data)
  rownames(testConsesus) <- rownames(data)
  countMat <- testConsesus;
  dataConcensus <- testConsesus
  totwts <- 0;
  for (i in 1:randomTests)
  {
    testset <- rownames(data[-randomSamples[[i]],])
    aclassLabels <- clusterLabels[[i]]$classification;
    nclus <- length(table(aclassLabels))
    wts <- (1.0-0.99*(nclus < 2))/(1.0+abs(nclus-numberofClusters));
    classLabels <- aclassLabels[testset];
    btestset <- rownames(data) %in% testset;
    for (id in testset)
    {
      testConsesus[id,btestset] <- testConsesus[id,btestset] + wts*(classLabels == aclassLabels[id]);
      countMat[id,btestset] <- countMat[id,btestset] + wts;
    }
    classLabels <- clusterLabels[[i]]$classification;
    for (id in 1:nrow(data))
    {
      dataConcensus[id,] <- dataConcensus[id,] + wts*(classLabels == classLabels[id]);
    }
    totwts <- totwts + wts;
  }
  cat("After Counting.")
  testConsesus[countMat > 0] <- testConsesus[countMat > 0]/countMat[countMat > 0];
  dataConcensus <- dataConcensus/totwts;
  pac <- sum(testConsesus[(testConsesus > pac.thr) & (testConsesus < (1.0 - pac.thr))])/nrow(data)/nrow(data);


  result <- list(randIndex = randIndex,jaccIndex = jaccIndex,meanJaccard = meanJaccard,randomSamples = randomSamples,
                 clusterLabels=clusterLabels, jaccardpoint=jaccardpoint,averageNumberofClusters=numberofClusters,
                 testConsesus=testConsesus,trainRandIndex = trainrandIndex,trainJaccIndex = trainjaccIndex,trainMeanJaccard = trainmeanJaccard,
                 trainJaccardpoint=trainjaccardpoint,PAC=pac,dataConcensus=dataConcensus);
  class(result) <- "ClusterStability"
  return(result);
}


clusterOutcomeStability <- function(clusterLabels=NULL,randomSamples=NULL,outcomeLabels=NULL)
{
  outcomRandIndex <- numeric();
  outcomejaccIndex <- numeric();
  outcomemeanJaccard <- numeric();
  for (i in 1:length(clusterLabels))
  {
    if (!is.null(outcomeLabels))
    {
      outcomRandIndex <- c(outcomRandIndex,adjustedRandIndex(clusterLabels[[i]]$classification[-randomSamples[[i]]],outcomeLabels[-randomSamples[[i]]]));
      outcomeJaccard <- jaccardMatrix(clusterLabels[[i]]$classification[-randomSamples[[i]]],outcomeLabels[-randomSamples[[i]]]);
      outcomejaccIndex <- c(outcomejaccIndex,outcomeJaccard$balancedMeanJaccard);
      outcomemeanJaccard <- c(outcomemeanJaccard,mean(outcomeJaccard$elementJaccard));
    }
  }

  result <- list(outcomRandIndex=outcomRandIndex,outcomejaccIndex=outcomejaccIndex,outcomemeanJaccard=outcomemeanJaccard);
  class(result) <- "ClusterOutcomeStability"
  return(result);
}

# Label the subjects that shere the same connectivity
getConsensusCluster <- function(object,who="training",thr=seq(0.80,0.30,-0.1))
{

  orgnames <-  rownames(object$dataConcensus);
  if (who != "training")
  {
    orgnames <-  rownames(object$testConsesus);
    pointJaccard <- object$jaccardpoint;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$testConsesus[order(-pointJaccard),]
  }
  else
  {
    pointJaccard <- 0.9*object$trainJaccardpoint + 0.1*object$jaccardpoint;
    names(pointJaccard) <- orgnames;
    concensusMat <- object$dataConcensus[order(-pointJaccard),]
  }
  concensusMat <- concensusMat[,order(-pointJaccard)]
  classID <- numeric(nrow(concensusMat));
  pointJaccard <- pointJaccard[order(-pointJaccard)];
  names(classID) <-  names(pointJaccard);
  npoints <- length(pointJaccard)
  label <- 1;
  for (lthr in thr)
  {
    totlabeled <- sum(classID > 0);
    if (totlabeled < npoints)
    {
      added <- 1;
      while (added > 0)
      {
        added <- 0;
        for (i in 1:npoints)
        {
          if (classID[i] == 0)
          {
              maxConnectedLabel <- label;
              wcon <- concensusMat[i,];
              consensA <- (wcon > lthr) & (classID > 0)
              consensB <- (wcon >= lthr) & (classID == 0)
              SconA <- sum(pointJaccard[consensA]*wcon[consensA]);
              SconB <- sum(pointJaccard[consensB]*wcon[consensB]) - pointJaccard[i]*wcon[i];
#              cat("A:",SconA,"B:",SconB,"P:",pointJaccard[i],"\n")

              if ( (SconB > 0.0075*npoints) || (SconA > 0.01*npoints) )
              {
                if (SconB > 0.75*SconA)
                {
                  classID[consensB] <- label;
                  added <- 1;
                  label <- label + 1;
                }
                else
                {
                  if (SconB >= 0)
                  {
                    if (SconA > 0)
                    {
                      tb <- table(classID[consensA])
                      smp <- tb;
                      for (nt in names(tb))
                      {
                        ptss <- consensA & (classID == as.numeric(nt));
                        smp[nt] <- sum(pointJaccard[ptss]*wcon[ptss]);
                      }
                      maxConnectedLabel <- names(which.max(smp))[1]
                      if ((0.75*smp[maxConnectedLabel]) <= SconB )
                      {
                        maxConnectedLabel <- label;
                      }
                      maxConnectedLabel <- as.numeric(maxConnectedLabel);
#                      print(c(label,smp,SconB,maxConnectedLabel))
                    }
                    classID[consensB] <- maxConnectedLabel;
                    added <- 1;
                    if (maxConnectedLabel == label)
                    {
                      label <- label + 1;
                    }
                  }
                }
              }
          }
        }
      }
      totlabeled <- sum(classID > 0);
      cat(maxConnectedLabel,":",sprintf("(%5.3f)",lthr),": ",totlabeled,": ")
    }
  }
  classID <- classID[orgnames];
  return (classID);
}

plot.ClusterStability <- function(object,...)
{
  plot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
  boxplot(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}

summary.ClusterStability <- function(object,...)
{
  summary(as.data.frame(cbind(randindex=object$randIndex,jaccardIndex=object$jaccIndex,meanJaccard=object$meanJaccard)),...)
}
