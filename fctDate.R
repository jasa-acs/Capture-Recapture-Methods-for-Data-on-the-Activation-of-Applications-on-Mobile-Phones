#' fctDate
#' 
#' Fonction qui ajoute un vecteur date  designant la date de chacune 
#' des occasions de capture de la matrice passee en parametre.  
#' 
#' 
#' @param matrice Dans le cas ou le parametre des=TRUE il sagit de la matrice 
#' de sortie de  sortie de la fonction \code{\link{descriptive}} du package 
#' \pkg{Rcapture} sinon il s'agit de la sortie de la fonction 
#' \code{\link{jollySeber}}. 
#' 
#' @param des Un objet de type logique. Valant TRUE lorsque matrice
#'  est la sortie de la fonction \code{\link{descriptive}} du package \pkg{Rcapture}  et
#'  FALSE lorsque matrice est la sortie de la fonction \code{\link{jollySeber}}.
#' 
#' @param jour Un objet de type logique. Valant TRUE lorsque les occasions 
#' de captures sont quotidiennes et FALSE lorsquelles sont hebdomadaire.  
#' 
#' @param date_debut Un objet de type format de date designant 
#' la date de debut du jeu de donnees.  
#' 
#' @return data_matrice Un objet de type data frame contenant 
#' la matrice de depart avec une colonne supplementaire designant la
#' date associee a chacune des occasions de capture. Dans le cas ou le parametre
#' des vaux FALSE les dates des occasions de capture suivent les indices des trois 
#' premiers vecteurs colonnes de la matrice (p, N et erreurType(N)). 
#' C'est-a-dire que le vecteur date dÃ©bute avec la date de la deuxieme 
#' occasion de capture et il termine avec la date de l'avant-derniere
#' occasion de capture. 
#' 
#' 
#' @author Joannie Houle
#'
fctDate <- function(matrice=matrice,date_debut=date_debut,des=FALSE, jour=TRUE)
{
  
  date1<-date_debut- 1
  if(jour==TRUE& des==FALSE)
  {
    n<-length(matrice[,1])
    vect_date<-date1 + seq(2, n+1, by = 1)
    
  } else if(jour==FALSE & des==FALSE) {
    n<-length(matrice[,1])
    vect_date<-date1 + seq(8, (n+1)*7, by = 7)
    
  } else if(jour==TRUE & des==TRUE){
    matrice<-matrice$base.freq
    n<-length(matrice[,1])
    vect_date<-date1 + seq(1, n, by = 1)
  } else {
    matrice<-matrice$base.freq
    n<-length(matrice[,1])
    vect_date<-date1 + seq(1, (n)*7, by = 7)
  }
  
  data_matrice<-as.data.frame(matrice)
  data_matrice$vect_date<-vect_date
  return(data_matrice)
  
}
