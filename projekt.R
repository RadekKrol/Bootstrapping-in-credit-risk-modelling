################################################################################
###                             Biblioteki                                   ###
################################################################################


rm(list = ls())

if(!require(R6)){install.packages("R6")  
  library(R6)}

if(!require(readr)){install.packages("readr")  
  library(readr)}

if(!require(tidyverse)){install.packages("tidyverse")  
  library(tidyverse)}

if(!require(caret)){install.packages("caret")  
  library(caret)}

if(!require(MASS)){install.packages("MASS")  
  library(MASS)}

if(!require(ModelMetrics)){install.packages("ModelMetrics")  
  library(ModelMetrics)}

if(!require(assertive)){install.packages("assertive") 
  library(assertive)
}
if(!require(ggplot2)){install.packages("ggplot2") 
  library(ggplot2)
}
if(!require(ROCit)){install.packages("ROCit") 
  library(ROCit)
}
if(!require(ggthemes)){install.packages("ggthemes") 
  library(ggthemes)
}


################################################################################
###                               Dane                                       ###
################################################################################

# podmienic sciezke
setwd("C:/Users/radkr/Desktop/bootstraping")
source(file = "pakiety.R")
source(file = "logit_bootstrap.R")
dane <- readr::read_csv2("german_credit_data.csv")
response <- readr::read_csv2("german_credit_response.csv")

################################################################################
###                             Obroka Danych                                ###
################################################################################

# Merge of 2 dataset
x <- merge.data.frame(dane, response, by.x = c("Age", "Credit amount", "Duration"), 
                      by.y = c("Age (years)", "Credit Amount", "Duration of Credit (month)"),
                      all.x = TRUE, all.y = FALSE )

# renaming and selecting columns
german <- x %>% dplyr::select("Age","Credit amount","Duration","Sex","Job","Housing",
                       "Saving accounts","Checking account", "Purpose.x", "Creditability" )

colnames(german) <- c("Age","Credit_amount","Duration","Sex","Job","Housing","Saving_acc",
                      "Checking_acc", "Purpose", "target" )

# String as factors
german[sapply(german, is.character)] <- lapply(german[sapply(german, is.character)], 
                                       as.factor)
german$Job <- as.factor(german$Job)

# Missing check
german <- german[!is.na(german$target),]

# For the rest of missing we create a dummy category
sum(is.na(german$target))

levels(german$Saving_acc) <- c("little","moderate","quite rich", "rich", "brak" )
levels(german$Checking_acc) <- c("little","moderate", "rich", "brak" )
german$Saving_acc[is.na(german$Saving_acc)]<- "brak"
german$Checking_acc[is.na(german$Checking_acc)]<- "brak"

length(german$target)
lapply(german, function(x){sum(is.na(x))})
glimpse(german)

# 1 oznacza default - kod do optymalizacji
german$target[german$target == 1] <- 2   
german$target[german$target == 0] <- 1   
german$target[german$target == 2] <- 0   

################################################################################
###                             Wstepny Model                                ###
################################################################################


# train- test split
set.seed(90)
in.train <- caret::createDataPartition(as.factor(german$target), p=0.6, list=FALSE)


train <- german[in.train,]
test <- german[-in.train,]


model <- glm(target~Age+Duration+Sex+Housing+Saving_acc+Checking_acc , family = "binomial", train)


gini(train$target, predict(model, train, type = "response"))


dum$glm
model

################################################################################
###                             Bootstrap                                    ###
################################################################################



logistic.model <- glm(target~Age+Duration+Sex+Housing+Saving_acc+Checking_acc , family = "binomial", train)

class <- logistic.model$y
score <- logistic.model$fitted.values

measure <- ROCit::measureit(score = score, class = class,
                     measure = c("ACC", "SENS", "FSCR"))



# accuracy plot

xxx$Cutoff[-1]
xxx$ACC[-1]
xxx <- as.data.frame(cbind(measure$ACC[-1], measure$Cutoff[-1]))
names(xxx) <- c("ACC", "Cutoff")
ggplot(xxx, aes(Cutoff, ACC)) + geom_line()+
  ggthemes::theme_pander(gM = TRUE, gm = FALSE, gl = "dotted", gc = "black" ) +
  xlab("Cut off level")+
  ylab("Accuracy")+
  ylim(c(0,1)) +
  geom_point(aes(x=Cutoff[which.max(ACC)], y=max(ACC)), colour="red", cex = 2) +
  annotate("text", x = xxx$Cutoff[which.max(xxx$ACC)], y = max(xxx$ACC)+ 0.05, 
           label = paste0("Best accuracy: ",round(max(xxx$ACC),2), " at cut-off: ", round(xxx$Cutoff[which.max(xxx$ACC)],2 )), 
           col = "red")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin=unit(c(2,1,2,1),"cm"))



# -----------------------------------------------
# ROC
# -----------------------------------------------






roc_empirical <- ROCit::rocit(score = score, class = class, negref = 0) 
AUC <- roc_empirical$AUC
gini <- 2*AUC - 1

### Roc curve & KS

y1 <- roc_empirical$TPR
x1 <- roc_empirical$FPR

index <- which.max(y1 - x1)

max(y1 - x1) # KS stat

roc_empirical$Cutoff[index]

plot(roc_empirical, values = F)
graphics::text(x = (x1[index] + 0.38), y = (y1[index] - 0.1),  as.character(round(roc_empirical$Cutoff[index],3)), 
               font = 3, cex = 0.95)




### Confidence interval plot - ROC curve

ciROC_emp90 <- ciROC(roc_empirical, 
                     level = 0.9)


col <- rep(c("#2F4F4F", "#404040"), 2)
lty <- rep(c(1, 2), 2)
lwd <- rep(c(2, 1), 2)
graphics::plot(ciROC_emp90$TPR ~ ciROC_emp90$FPR, xlab = "1 - Specificity (FPR)", 
               ylab = "Sensitivity (TPR)", type = "l", col = col[1], 
               lwd = lwd[1], lty = lty[1])
graphics::grid()
graphics::lines(ciROC_emp90$LowerTPR ~ ciROC_emp90$FPR, lty = lty[2], lwd = lwd[2], 
                col = col[2])
graphics::lines(ciROC_emp90$UpperTPR ~ ciROC_emp90$FPR, lty = lty[2], lwd = lwd[2], 
                col = col[2])
legend("bottomright", c("ROC curve","90% Confidence Interval"),
       lty = c(1,2), col = c(1,1), lwd = c(2,1))


### KS plot 


kplot <- ROCit::ksplot(roc_empirical)

### score distributions

temp <- as.data.frame(logistic.model$data$target)
colnames(temp) <- "target"

temp["score"] <- predict(logistic.model, logistic.model$data, type = "response")
head(temp)


temp %>% 
  group_by(target) %>% 
  summarise(tb = mean(score)) %>% 
  ungroup() -> mean_score_train


temp %>% 
  ggplot(aes(x=score, fill = factor(target))) + 
  geom_density(alpha = 0.3, aes(group = target))+
  geom_vline(aes(xintercept = mean_score_train$tb[1]), linetype = "dashed", color = "red") + 
  geom_vline(aes(xintercept = mean_score_train$tb[2]), linetype = "dashed", color = "blue") +
  geom_text(aes(x = mean_score_train$tb[1]+ 0.01, y = max(density(temp$score)$y), label = mean_score_train$tb[1] %>% round(2)), color = "red", size = 4) + 
  geom_text(aes(x = mean_score_train$tb[2]+ 0.01, y = max(density(temp$score)$y), label = mean_score_train$tb[2] %>% round(2)), color = "blue", size = 4) +
  theme_classic()+
  theme(legend.position = "top" ) +
  scale_fill_discrete(name = "model_score", labels = c(0, 1))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin=unit(c(2,1,2,1),"cm"))




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

### R6 class
logit_bootstrap <- R6::R6Class(classname = "logit_Bootstrap",
                               
  # ----------------------------------------------------------------------------
  # private
  # ----------------------------------------------------------------------------          
  
  private = list(
    ..seed = 11,
    ..model = NULL,
    ..data = NULL,
    ..target = NULL,
    ..train = NULL,
    ..test = NULL,
    ..mat_coef = NULL,
    ..results = NULL,
    ..is_data_split = FALSE,
    ..boot_model = NULL,
    ..stat = NULL
  ),
  
  # ----------------------------------------------------------------------------
  # public
  # ----------------------------------------------------------------------------  
  
  public = list(
    
    
    split = function(p = 0.6){
      
      if (!(requireNamespace("caret"))){
        stop("The 'caret' package is missing")
      }
      
      if(p>=1 || p <= 0){
        stop("Argument p musi zawiera� si� w przedziale (0,1)")
      }
      
      `%G%` <- paste0
      set.seed(private$..seed)
      in.train <- caret::createDataPartition(as.factor(private$..data[,private$..target]), p=p, list=FALSE)
      private$..train <- private$..data[in.train,]
      private$..test <- private$..data[-in.train,]
      message("Dokonano podzia�u danych na zbi�r treningowy i walidacyjny w stosunku " %G% as.character(p*10) %G% ":" %G% as.character(10-p*10))
      private$..is_data_split <- TRUE
      
    },
    
    unsplit = function(){
      
      message("Podzia� na zbi�r trenngowy i walidacyjny zosta� usuni�ty")
      private$..train <- NULL
      private$..test <- NULL
      private$..is_data_split <- FALSE
      
    },
    
    bootstrap = function(n_boot = 100){
      
      
      if((!is.numeric(n_boot)) || n_boot < 1 || n_boot > 1000000){
        stop("Argument n_boot powinien przyjmowa� liczby ca�kowite z zakresu <1,1 000 000>")
      }
      
      n <- length(names(private$..model$coefficients))
      
      private$..mat_coef <- matrix(NA, nrow=n_boot, ncol = n)
      
      colnames(private$..mat_coef) <- names(private$..model$coefficients)
      
      private$..stat <- matrix(NA, nrow=n_boot, ncol = 2)

      if(private$..is_data_split == FALSE){
        
        message("Przeprowadzanie bootstrapingu na pe�nych danych")
        
        for(i in 1:n_boot){
          dane <- sample(length(private$..data[,1]), replace = TRUE)
          m1 <- glm(formula = private$..model$formula, family = "binomial", data=private$..data[dane,])
          private$..mat_coef[i,] <- m1$coefficients
          roc_empirical <- ROCit::rocit(score = m1$fitted.values, class = m1$y, negref = 0) 
          private$..stat[i,1] <- 2*roc_empirical$AUC - 1
          private$..stat[i,2] <- max(roc_empirical$TPR - roc_empirical$FPR)
        }
      } else {
        
        message("Przeprowadzanie bootstrapingu na zbiorze treningowym")
        
        for(i in 1:n_boot){
          dane <- sample(length(private$..train[,1]), replace = TRUE)
          m1 <- glm(formula = private$..model$formula, family = "binomial", data=private$..train[dane,])
          private$..mat_coef[i,] <- m1$coefficients
          roc_empirical <- ROCit::rocit(score = m1$fitted.values, class = m1$y, negref = 0) 
          private$..stat[i,1] <- 2*roc_empirical$AUC - 1
          private$..stat[i,2] <- max((roc_empirical$TPR - roc_empirical$FPR))
        }
      }
      
      stats <- matrix(NA, nrow = n, ncol = 16)
      rownames(stats) <- colnames(private$..mat_coef)
      colnames(stats) <- c("original", "bootEst", "bootBias", 
                           "bootSE", "bootMed", "bootMin", "bootMax", 
                           "1%",  "5%",  "10%", "25%", "50%", "75%", "90%", "95%", "99%")
      
      stats[, 1] <- private$..model$coefficients
      stats[, 2] <- apply(private$..mat_coef, 2, function(x) mean(x, na.rm = TRUE))
      stats[, 3] <- stats[, 2] - stats[, 1]
      stats[, 4] <- apply(private$..mat_coef, 2, function(x) median(x, na.rm = TRUE))
      stats[, 5] <- apply(private$..mat_coef, 2, function(x) sd(x, na.rm = TRUE))
      stats[, 6] <- apply(private$..mat_coef, 2, function(x) min(x, na.rm = TRUE))
      stats[, 7] <- apply(private$..mat_coef, 2, function(x) max(x, na.rm = TRUE))
      stats[, 8:16] <- t(apply(private$..mat_coef, 2, 
                               function(x) quantile(x, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)))

      
      private$..stat <- as.data.frame(private$..stat)
      private$..results <- as.data.frame(stats)
      
      
      # oszacowanie ko�cowego modelu ze �rednimi wsp�czynnikami
      private$..boot_model <- private$..model
      private$..boot_model$coefficients <- apply(private$..mat_coef, 2, mean)
      names(private$..boot_model$coefficients) <- names(private$..model$coefficients)
      
      
      return(private$..results[1:7])
    },
    
    # what = c("coef", "gini", "KS")
    conf_intervals = function(what = "coef", simplify = TRUE){
      
      if(!is.logical(simplify)){
        stop("Argument simplify musi mie� typ logiczny")
      }

      
      
      if(is.null(private$..results)){
        stop("Najpierw nale�y oszacowa� bootstrapowy model")
      }
      if(!(what %in% c("coef", "gini", "KS"))){
        stop("Niepropwana warto�� argumentu WHAT = c('coef', 'gini', 'KS')")
      }
      
      
      if(what == "coef"){
        
        message("Confidence Interval - Coefficients")
        
        if(simplify){
          return(private$..results[,c(9,15)])
        } else {
          return(private$..results[,8:16])
        }
      } else if(what == "gini"){
        
        message("Confidence Interval - Gini")
        message("Statystyk� oszacowane na zbiorze treningowym")
        if(simplify){
          return(quantile(private$..stat[,1], probs = c(0.05, 0.95), na.rm = TRUE))
        } else {
          return(quantile(private$..stat[,1], probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE))
        }
      } else if(what == "KS"){
        
        message("Confidence Interval - KS")
        message("Statystyk� oszacowane na zbiorze treningowym")
        
        if(simplify){
          return(quantile(private$..stat[,2], probs = c(0.05, 0.95), na.rm = TRUE))
        } else {
          return(quantile(private$..stat[,2], probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE))
        }
      }
    },
    
    discriminatory_power = function(boot = TRUE, test_data = TRUE){
      
      if(!is.logical(boot)){
        stop("Argument boot musi mie� typ logiczny")
      }
      if(!is.logical(test_data)){
        stop("Argument test_data musi mie� typ logiczny")
      }
      
      
      if(boot){
        if(is.null(private$..boot_model)){
          stop("Nie przeprowadzono dotychczas procedury bootstrapu")
        }
        message("Statystyki dyskryminacyjne modelu bootstrapowego")
        model <- private$..boot_model
      } else {
        message("Statystyki dyskryminacyjne modelu pocz�tkowego")
        model <- private$..model
      }
      
      if(test_data){
        if(!private$..is_data_split){
          stop("Nie wydzielono zbioru walidacyjnego")
        }
        message("Statystyki dyskryminacyjne modelu wyliczone na zbiorze walidacyjnym")
        dt <- private$..test
      } else {
        message("Statystyki dyskryminacyjne modelu wyliczone na zbiorze treningowym")
        if(!private$..is_data_split){
          dt <- private$..data
        }else{
        dt <- private$..train
        }
      }
      
      
      mylist <- sapply(c("gini", "AUC", "KS", "Optimal cutoff - Younden", "Optimal cutoff - accuracy",
                         "metrics for optimal cutoff - Younden", "metrics for optimal cutoff - accuracy"),
                       function(x) NULL)
      
      roc_empirical <- ROCit::rocit(score = predict(model, dt, type = "response"), class = dt[,private$..target], negref = 0) 
      index <- which.max((roc_empirical$TPR - roc_empirical$FPR))
      
      measure <- ROCit::measureit(score = predict(model, dt, type = "response"), class = dt[,private$..target],
                                  measure = c("ACC"))
      

      mylist[[1]] <- gini(dt[,private$..target], predict(model, dt, type = "response"))
      mylist[[2]] <- (mylist[[1]] + 1)/2
      mylist[[3]] <- max((roc_empirical$TPR - roc_empirical$FPR))
      mylist[[4]] <- roc_empirical$Cutoff[index]
      mylist[[5]] <- measure$Cutoff[which.max(measure$ACC)]
      
      pred_Y <- as.factor(ifelse(predict(model, dt, type="response")>roc_empirical$Cutoff[index],"1","0"))
      pred_A <- as.factor(ifelse(predict(model, dt, type="response")>measure$Cutoff[which.max(measure$ACC)],"1","0"))
      actual <- as.factor(dt[,private$..target])
        
      mylist[[6]] <- caret::confusionMatrix(pred_Y, actual)
      mylist[[7]] <- caret::confusionMatrix(pred_A, actual)
      
      
      return(mylist)
      
    },
    
    # what = c("accuracy", "ROC", "ROC_CI", "KS", "Score")
    # level wyznacza przedzia� ufno�ci je�li rysujesz ROC_CI
    make_plots = function(what = "accuracy", boot = TRUE, test_data = TRUE, level = 0.9){
      
      if(level >= 1 || level <= 0){
        stop("Argument level musi zawiera� si� w przedziale (0,1)")
      }
      
      if(!(what %in% c("accuracy", "ROC", "ROC_CI", "KS", "Score"))){
        stop("Niepropwana warto�� argumentu WHAT = c('coef', 'gini', 'KS')")
      }
      
      if(!is.logical(boot)){
        stop("Argument boot musi mie� typ logiczny")
      }
      if(!is.logical(test_data)){
        stop("Argument test_data musi mie� typ logiczny")
      }
      
      if(boot){
        if(is.null(private$..boot_model)){
          stop("Nie przeprowadzono dotychczas procedury bootstrapu")
        }
        message("Wizaulizacja dla modelu bootstrapowego")
        model <- private$..boot_model
      } else {
        message("Wizaulizacja dla modelu modelu pocz�tkowego")
        model <- private$..model
      }
      
      if(test_data){
        if(!private$..is_data_split){
          stop("Nie wydzielono zbioru walidacyjnego")
        }
        message("Wizaulizacja dla zbioru walidacyjnego")
        dt <- private$..test
      } else {
        message("Wizaulizacja dla zbioru treningowego")
        if(!private$..is_data_split){
          dt <- private$..data
        }else{
          dt <- private$..train
        }
      }
      
      
      roc_empirical <- ROCit::rocit(score = predict(model, dt, type = "response"), 
                                    class = dt[,private$..target], negref = 0) 
      
      
      if(what == "accuracy"){
        
        measure <- ROCit::measureit(score = predict(model, dt, type = "response"), 
                                    class = dt[,private$..target],
                                    measure = c("ACC"))

        xxx <- as.data.frame(cbind(measure$ACC[-1], measure$Cutoff[-1]))
        names(xxx) <- c("ACC", "Cutoff")
        ggplot(xxx, aes(Cutoff, ACC)) + geom_line()+
          ggthemes::theme_pander(gM = TRUE, gm = FALSE, gl = "dotted", gc = "black" ) +
          xlab("Cut off level")+
          ylab("Accuracy")+
          ylim(c(0,1)) +
          geom_point(aes(x=Cutoff[which.max(ACC)], y=max(ACC)), colour="red", cex = 2) +
          annotate("text", x = xxx$Cutoff[which.max(xxx$ACC)], y = max(xxx$ACC)+ 0.05, 
                   label = paste0("Best accuracy: ",round(max(xxx$ACC),2), " at cut-off: ", round(xxx$Cutoff[which.max(xxx$ACC)],2 )), 
                   col = "red")+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                plot.margin=unit(c(2,1,2,1),"cm"))
        
        
      }else if(what == "ROC"){
        
        index <- which.max((roc_empirical$TPR - roc_empirical$FPR))
        
        plot(roc_empirical, values = F)
        graphics::text(x = roc_empirical$FPR[index] + 0.38, y = roc_empirical$TPR[index] - 0.1,  
                       as.character(round(roc_empirical$Cutoff[index],3)), 
                       font = 3, cex = 0.95)
        
        
        
      }else if(what == "ROC_CI"){
        
        ciROC_emp90 <- ciROC(roc_empirical, 
                             level = level)
        
        
        col <- rep(c("#2F4F4F", "#404040"), 2)
        lty <- rep(c(1, 2), 2)
        lwd <- rep(c(2, 1), 2)
        graphics::plot(ciROC_emp90$TPR ~ ciROC_emp90$FPR, xlab = "1 - Specificity (FPR)", 
                       ylab = "Sensitivity (TPR)", type = "l", col = col[1], 
                       lwd = lwd[1], lty = lty[1])
        graphics::grid()
        graphics::lines(ciROC_emp90$LowerTPR ~ ciROC_emp90$FPR, lty = lty[2], lwd = lwd[2], 
                        col = col[2])
        graphics::lines(ciROC_emp90$UpperTPR ~ ciROC_emp90$FPR, lty = lty[2], lwd = lwd[2], 
                        col = col[2])
        legend("bottomright", c("ROC curve",paste(round(level*100),"% Confidence Interval")),
               lty = c(1,2), col = c(1,1), lwd = c(2,1))
        
      }else if(what == "KS"){
        
        ROCit::ksplot(roc_empirical)
        
      }else if(what == "Score"){
        
        temp <- as.data.frame(dt[,private$..target])
        colnames(temp) <- "target"
        
        temp["score"] <- predict(model, dt, type = "response")
        
        
        temp %>% 
          group_by(target) %>% 
          summarise(tb = mean(score)) %>% 
          ungroup() -> mean_score_train
        
        
        temp %>% 
          ggplot(aes(x=score, fill = factor(target))) + 
          geom_density(alpha = 0.3, aes(group = target))+
          geom_vline(aes(xintercept = mean_score_train$tb[1]), linetype = "dashed", color = "red") + 
          geom_vline(aes(xintercept = mean_score_train$tb[2]), linetype = "dashed", color = "blue") +
          geom_text(aes(x = mean_score_train$tb[1]+ 0.01, y = max(density(temp$score)$y), label = mean_score_train$tb[1] %>% round(2)), color = "red", size = 4) + 
          geom_text(aes(x = mean_score_train$tb[2]+ 0.01, y = max(density(temp$score)$y), label = mean_score_train$tb[2] %>% round(2)), color = "blue", size = 4) +
          theme_classic()+
          theme(legend.position = "top" ) +
          scale_fill_discrete(name = "model_score", labels = c(0, 1))+
          theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                plot.margin=unit(c(2,1,2,1),"cm"))
        
      }
    
    },
    
    # # what = c("coef", "gini", "KS")
    show_density = function(){
      
      if(is.null(private$..results)){
        stop("Najpierw nale�y oszacowa� bootstrapowy model")
      }
      if(!(what %in% c("coef", "gini", "KS"))){
        stop("Niepropwana warto�� argumentu WHAT = c('coef', 'gini', 'KS')")
      }
      
      
      g <- reshape2::melt(as.data.frame(private$..mat_coef))
      
      ggplot(data = g, aes(value))+geom_density()+
        geom_density(alpha = 0.3, fill="purple", color = "purple")+
        facet_wrap(~variable, scales = "free")+
        geom_vline(aes(xintercept = 0), linetype = "solid", color = "#FF6666", size =1) +
        theme_classic()
      
    },
    
    initialize = function(model, data, target, seed){
      
      
      # Sprawdzenie czy przekazano dany parametr
      if(!missing(seed)){ 
        private$..seed <- seed
      }
      if(missing(target)){ 
        stop("Nie okre�lono zmiennej celu")
      } else if (missing(data)){
        stop("Nie wskazano zbioru danych")
      } else if(missing(model)) {
        stop("Nie wskazano obiektu glm")
      }
      
      # czy dane s� obiektem klas data.frame
      if(!assertive::is_data.frame(data)){
        stop("Przekazany zbi�r danych musi by� obiektem typu data.frame")
      } else {
        private$..data <- as.data.frame(data)
      }
      
      # czy zmienna celu znajduje si� w danych
      if(!assertive::is_character(target)){
        stop("Zmienna celu musi mie� typ character")
      } else if(!(target %in% names(private$..data))) {
        stop("Brak zmiennej celu w przekazanym zbiorze danych")
      } else {  
        private$..target <- target
      }
      
      # Czy przekazano poprawny model
      objasniana <- model$formula[2][[1]]
      objasniajace <- unlist(strsplit(deparse(model$formula[3][[1]]), " \\+ "))
      
      check1 <- !(deparse(objasniana) %in% names(data))
      
      check2 <- !(any(unlist((lapply(objasniajace, 
                         function(x){x %in% names(data)})))))
      
      if(!(class(model)[1] == "glm")){
        stop("Nale�y przekaza� obiekt o klasie 'glm'")
      } else if(check1 || check2){
        stop("W danych nie odnaleziono wszystkich zmiennych zadanych w modelu")
      } else {
        private$..model <- model
      }
      
    }
  ),

  # ----------------------------------------------------------------------------
  # active
  # ----------------------------------------------------------------------------  
  
  active = list(
    
    # seed = function(x){
    #   if(missing(x)){
    #     return(private$..seed)
    #   } else if(!is.numeric(x) || (trunc(x) != x) ){
    #     message("Nie uda�o si� zmieni� seeda - nale�y poda� liczb� ca�kowit�")
    #   } else {  
    #     private$..seed <- x
    #   }
    # },
    
    data = function(x){
      if(missing(x)){
        return(private$..data)
      } else {
        stop("Nie mo�na nadpisa� tego atrybutu")
      }
    },
    
    model = function(x){
      if(missing(x)){
        return(private$..model$formula)
      } else {
        stop("Nie mo�na nadpisa� tego atrybutu")
      }
    },
    
    results = function(x){
      if(missing(x)){
        return(private$..results)
      } else {
        stop("Nie mo�na nadpisa� tego atrybutu")
      }
    },
    
    is_data_split = function(x){
      if(missing(x)){
        return(private$..is_data_split)
      } else {
        stop("Nie mo�na nadpisa� tego atrybutu")
      }
    }
    
  )
)



dum <- logit_bootstrap$new(model, as_tibble(german), "target", seed = 1000001)

dum$split()

nn<- dum$bootstrap(n_boot = 100)

dum$make_plots(what = "Score")

dum$show_distribution()

a <- dum$discriminatory_power(boot = TRUE, test_data = TRUE)
b<- dum$discriminatory_power(boot = FALSE, test_data = TRUE)
c<- dum$discriminatory_power(boot = TRUE, test_data = FALSE)
d<-dum$discriminatory_power(boot = FALSE, test_data = FALSE)

a$gini
b$gini
c$gini
d$gini





########

# 1. wstep - opis problemu - zagadnienia bankowe

# 2. Opis zbioru danych

# 3. Opis dzialania bootstrapu

# 4. Przeksztalcenia danych (krotki fragment)

# 5. Opis wstepnego modelu oraz selekcja zmiennych


# 6.1 model obiektowy

# 6. Zastosowanie bootstrapu

# 7. Opisanie parametrow, bootstrapowych przedzialow ufnosci, oraz bootstrapowych statystyk
#    Porownanie ich z modelem pierwotnym oraz wizualizacja wynikow

# 8. Zakonczenie i podsumowanie 




