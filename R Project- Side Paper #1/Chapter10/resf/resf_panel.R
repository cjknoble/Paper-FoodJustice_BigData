resf_panel 	<- function( y, x = NULL, t_id, s_id, meig, pmodel = "random", effect = "twoways" ){
    
    if( pmodel == "pooling" ){
    		meig$sf<- meig$sf[ meig$other$s_id, ][ s_id, ]
    		meig$other$fast	<- 1
    		res		<- resf_p( y = y, x = x, meig = meig )
    		b_par	<- res$b
    		bg_par	<- NULL
    		s_par	<- res$s
    		sg_par	<- NULL
    		e_stat	<- res$e
    		r_par	<- res$r
    		SF		<- res$sf
    		pred	<- res$pred
    		resid	<- res$resid
    		other	<- res$other 

    } else if( pmodel == "within" ){
      x     <- data.frame( x )
      nx    <- dim( x )[2]
      if( effect == "individual" ){
    		s_id0 <- unique( s_id )
    		s_cent<- as.data.frame(matrix(0, nrow = length( s_id0 ), ncol = nx + 2)) # you must use the as.data.frame() routine to change the matrix of s_cent to a dataframe so that the s_id0 can be read as characters (if they are numeric ids, then it will be fine). I have to change this to all three _cent: the t_cent and the twoways s_cent. Consequentially, when calculating s_ef, s_cent[,"y"] has to be changed back to the matrix form using as.matrix.
    		s_cent[ ,1 ] <- s_id0
    		s_cent[ ,2 ] <- tapply(y, s_id, mean)
    		for( i  in 1:nx ) s_cent[ ,2 + i ]<- tapply( x[ ,i ], s_id, mean )
    		t_cent<- NULL
    		g_cent<- NULL
    		
    	  xg00	  <- factor( s_id )
    		xg0	    <- model.matrix( ~ 0 + xg00 )
    		xg0s	  <- colSums( xg0 )
    		xg_idid	<- max( which( xg0s == max( xg0s ) ) )
    		xg	    <- data.frame( xg0[ , -xg_idid ] )
    		names( xg ) <- paste( "s_", levels( xg00 )[ -xg_idid ], sep = "" )

      } else if( effect == "time" ){
        t_id0 <- unique( t_id )
        t_cent<- as.data.frame(matrix(0, nrow = length( t_id0 ), ncol = nx + 2))
        t_cent[ ,1 ] <- t_id0
        t_cent[ ,2 ] <- tapply(y, t_id, mean)
        for( i  in 1:nx ) t_cent[ ,2 + i ]<- tapply( x[ ,i ], t_id, mean )
        s_cent <- NULL
        g_cent <- NULL
        
        xg00	  <- factor( t_id )
        xg0	    <- model.matrix( ~ 0 + xg00 )
        xg0s	  <- colSums( xg0 )
        xg_idid	<- max( which( xg0s == max( xg0s ) ) )
        xg	    <- data.frame( xg0[ , -xg_idid ] )
        names( xg ) <- paste( "t_", levels( xg00 )[ -xg_idid ], sep = "" )
        
      } else if( effect == "twoways" ){
        s_id0 <- unique( s_id )
        s_cent<- as.data.frame(matrix(0, nrow = length( s_id0 ), ncol = nx + 2))
        s_cent[ ,1 ] <- s_id0
        s_cent[ ,2 ] <- tapply(y, s_id, mean)
        for( i  in 1:nx ) s_cent[ ,2 + i ]<- tapply( x[ ,i ], s_id, mean )
        
        t_id0 <- unique( t_id )
        t_cent<- as.data.frame(matrix(0, nrow = length( t_id0 ), ncol = nx + 2))
        t_cent[ ,1 ] <- t_id0
        t_cent[ ,2 ] <- tapply(y, t_id, mean)
        for( i  in 1:nx ) t_cent[ ,2 + i ]<- tapply( x[ ,i ], t_id, mean )
        
        xg00	  <- factor( s_id )
        xg0	    <- model.matrix( ~ 0 + xg00 )
        xg0s	  <- colSums( xg0 )
        xg_idid	<- max( which( xg0s == max( xg0s ) ) )
        xg_s	    <- data.frame( xg0[ , -xg_idid ] )
        names( xg_s ) <- paste( "s_", levels( xg00 )[ -xg_idid ], sep = "" )
        
        xg00	  <- factor( t_id )
        xg0	    <- model.matrix( ~ 0 + xg00 )
        xg0s	  <- colSums( xg0 )
        xg_idid	<- max( which( xg0s == max( xg0s ) ) )
        xg_t	  <- data.frame( xg0[ , -xg_idid ] )
        names( xg_t ) <- paste( "t_", levels( xg00 )[ -xg_idid ], sep = "" )
        xg      <- data.frame( xg_s, xg_t )
        
        g_cent <-c(mean(y), colMeans(x) )
      }
      
      y2    <- y
      x2    <- x
      g_np  <- 0
      if( is.null( s_cent ) == FALSE ){
        s_cent<-as.data.frame( s_cent )
        names( s_cent )[1] <- "s_id"
        names( s_cent )[2] <- "y"
        names( s_cent )[-(1:2)]<-names( x )
        
        s_dat  <- data.frame( s_id, id = 1:length( y ) )
        s_cent1<- merge( s_dat, s_cent, by = "s_id", all.x = TRUE)
        s_cent1<- s_cent1[ order( s_cent1$id ), ]
        y2     <- y2 - s_cent1$y
        x2     <- x2 - s_cent1[,-(1:3)]
        g_np   <- g_np + dim( s_cent )[1] - 1
      } else {
        s_cent1<- NULL
      }

      if( is.null( t_cent ) == FALSE ){
        t_cent<-as.data.frame( t_cent )
        names( t_cent )[1] <- "t_id"
        names( t_cent )[2] <- "y"
        names( t_cent )[-(1:2)]<-names( x )
        
        t_dat  <- data.frame( t_id, id = 1:length( y ) )
        t_cent1<- merge( t_dat, t_cent, by = "t_id", all.x = TRUE)
        t_cent1<- t_cent1[ order( t_cent1$id ), ]
        y2     <- y2 - t_cent1$y
        x2     <- x2 - t_cent1[,-(1:3)]
        g_np   <- g_np + dim( t_cent )[1] - 1
      } else {
        t_cent1<- NULL
      }

      if( is.null( g_cent ) == FALSE ){
        g_cent<-as.data.frame( g_cent )
        y2     <- y2 + g_cent[1,1]
        for( i in 1:length(g_cent[-1,1])  ){
          x2[ ,i ]     <- x2[ ,i ] + g_cent[ i + 1, 1 ]
        }
        g_np   <- g_np - 1
      }
      
      res		<- resf_p( y = y2, x = x2, y0 = y, meig = meig, s_id = s_id, g_np = g_np )
    	b_par	<- res$b
    	s_par	<- res$s
    	e_stat<- res$e
    	sg_par<- NULL
    	r_par	<- res$r
    	SF		<- res$sf
    	resid	<- res$resid
    	pred  <- y - resid
    	other	<- res$other
    	par00 <- other$par0
    	
    	const <- 0
    	bg_par<-NULL
    	if( is.null( s_cent1 ) == FALSE ){
    		  s_ef <- as.matrix(s_cent[, "y" ]) - as.matrix( s_cent [,-( 1:2 )]) %*% b_par[-1,1] # Here, s_cent[,"y"] must be changed to a matrix form before it can be used for calculating s_ef. This will be the same for the t_cent and s_cent in the twoways.
          s_mean<-mean( s_ef )
          s_ef <- s_ef - s_mean
          const<- const + s_mean
          
          bg_s<-data.frame(Estimate = s_ef)
          row.names(bg_s) <- paste( "s_", s_id0, sep="")
          bg_par<-rbind( bg_par, bg_s)
      }
    	if( is.null( t_cent1 ) == FALSE ){
    		  t_ef <- as.matrix(t_cent[ ,"y" ]) - as.matrix( t_cent [,-( 1:2 )]) %*% b_par[-1,1]
    		  t_mean<-mean( t_ef )
    		  t_ef <- t_ef - t_mean
    		  const<- const + t_mean

    		  bg_t<-data.frame(Estimate = t_ef)
    		  row.names(bg_t) <- paste( "t_", t_id0, sep="")
    		  bg_par <- rbind( bg_par, bg_t)
    	}
    		
    	if( is.null( g_cent ) == FALSE ){
    		  g_ef <- g_cent[ 1, 1 ] - c( g_cent [ -1, 1 ] %*% b_par[-1,1] )
    		  const<- const - g_ef
    	}
    	 b_par[ 1,1 ]<-const
    	 b_par[ 1, -1 ]<-NA

    	} else if( pmodel == "random" ){
    		if( effect == "individual" ){
    			xgroup	<- data.frame( s_id )
        } else if( effect == "time" ){
    			xgroup	<- data.frame( t_id )
    		} else if( effect == "twoways" ){
    		    xgroup	<- data.frame( s_id, t_id )
    		}
		    res		<- resf_p( y = y, x = x, meig = meig, xgroup = xgroup, s_id =s_id )
    		b_par	<- res$b
    		bg_par	<- res$b_g
    		s_par	<- res$s
    		sg_par	<- res$s_g
    		e_stat	<- res$e
    		r_par	<- res$r
    		SF		<- res$sf
    		pred	<- res$pred
    		resid	<- res$resid
    		other	<- res$other 
    }
    
    return( list( b = b_par, b_eff = bg_par, s = s_par, s_eff = sg_par, e = e_stat,
    		  r = r_par, sf = SF, pred = pred, resid = resid, other = other ) )
}
