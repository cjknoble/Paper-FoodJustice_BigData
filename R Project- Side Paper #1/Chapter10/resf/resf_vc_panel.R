resf_vc_panel		<- function( y, x = NULL, xconst = NULL, t_id, s_id, meig, 
			pmodel = "random", effect = "twoways", 
			penalty = "bic", maxiter = 30, sizelimit = 2000 ){

    if( pmodel == "pooling" ){
    		meig$sf<- meig$sf[ meig$other$s_id, ][ s_id, ]
    		meig$other$fast	<- 1
    		res		<- resf_vc_p( y = y, x = x, xconst = xconst, meig = meig, s_id = s_id,
    							penalty = penalty, maxiter = maxiter, sizelimit = sizelimit )
    		b_par	<- res$b
    		bg_par	<- NULL
    		s_par	<- res$s
    		sg_par	<- NULL
    		e_stat	<- res$e
    		b_vc		<- res$b_vc
    		bse_vc	<- res$bse_vc
    		t_vc		<- res$t_vc
    		p_vc		<- res$p_vc
    		pred	<- res$pred
    		resid	<- res$resid
    		vc		<- res$vc
    		r		<- res$r
    		other	<- res$other 

    } else if( pmodel == "within" ){
      x     <- data.frame( x )
      nx    <- dim( x )[2]
      if( is.null( xconst ) ==FALSE ) {
        xconst  <- data.frame( xconst )
        nxc     <- dim( xconst )[ 2 ]
      }
      sc_cent<- NULL
      tc_cent<- NULL
      gc_cent<- NULL
      
      if( effect == "individual" ){
        s_id0 <- unique( s_id )
        s_cent<- as.data.frame(matrix(0, nrow = length( s_id0 ), ncol = nx + 2)) # must use the as.data.frame routine to make sure that the next statement s_cent[,1] <-s_id0 is correctly receiving the names from the data. Because matrix will not hold non-numeric values. This will have to be changed for the time effect and twoways as well.
        s_cent[ ,1 ] <- s_id0
        s_cent[ ,2 ] <- tapply(y, s_id, mean)
        for( i  in 1:nx ) s_cent[ ,2 + i ]<- tapply( x[ ,i ], s_id, mean )
        t_cent<- NULL
        t_cent1<-NULL
        g_cent<- NULL

        if( is.null( xconst ) ==FALSE ){
          sc_cent<- matrix(0, nrow = length( s_id0 ), ncol = nxc + 2)
          sc_cent[ ,1 ] <- s_id0
          for( i  in 1:nxc ) sc_cent[ ,2 + i ]<- tapply( xconst[ ,i ], s_id, mean )
        }
                
      } else if( effect == "time" ){
        t_id0 <- unique( t_id )
        t_cent<- as.data.frame(matrix(0, nrow = length( t_id0 ), ncol = nx + 2))
        t_cent[ ,1 ] <- t_id0
        t_cent[ ,2 ] <- tapply(y, t_id, mean)
        for( i  in 1:nx ) t_cent[ ,2 + i ]<- tapply( x[ ,i ], t_id, mean )
        s_cent <- NULL
        s_cent1<- NULL
        g_cent <- NULL

        if( is.null( xconst ) ==FALSE ){
          tc_cent<- matrix(0, nrow = length( t_id0 ), ncol = nxc + 2)
          tc_cent[ ,1 ] <- t_id0
          for( i  in 1:nxc ) tc_cent[ ,2 + i ]<- tapply( xconst[ ,i ], t_id, mean )
        }
        
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
        
        if( is.null( xconst ) ==FALSE ){
          sc_cent<- as.data.frame(matrix(0, nrow = length( s_id0 ), ncol = nxc + 2))
          sc_cent[ ,1 ] <- s_id0
          for( i  in 1:nxc ) sc_cent[ ,2 + i ]<- tapply( xconst[ ,i ], s_id, mean )

          tc_cent<- as.data.frame(matrix(0, nrow = length( t_id0 ), ncol = nxc + 2))
          tc_cent[ ,1 ] <- t_id0
          for( i  in 1:nxc ) tc_cent[ ,2 + i ]<- tapply( xconst[ ,i ], t_id, mean )
          
          gc_cent<- data.frame(c(mean(y), colMeans(xconst)))
        }
        
        g_cent <-c(mean(y), colMeans(x) )
      }
      
      y2    <- y
      x2    <- x
      xconst2<-xconst
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
        
        if( is.null( xconst ) == FALSE){
          sc_cent<-as.data.frame( sc_cent )
          names( sc_cent )[1] <- "s_id"
          names( sc_cent )[-(1:2)]<-names( xconst )
          sc_cent1<- merge( s_dat, sc_cent, by = "s_id", all.x = TRUE)
          sc_cent1<- sc_cent1[ order( sc_cent1$id ), ]
          xconst2 <- xconst2 - sc_cent1[,-(1:3)]
        }        
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

        if( is.null( xconst ) == FALSE){
          tc_cent<-as.data.frame( tc_cent )
          names( tc_cent )[1] <- "t_id"
          names( tc_cent )[-(1:2)]<-names( xconst )
          tc_cent1<- merge( t_dat, tc_cent, by = "t_id", all.x = TRUE)
          tc_cent1<- tc_cent1[ order( tc_cent1$id ), ]
          xconst2 <- xconst2 - tc_cent1[,-(1:3)]
        }
      }
      
      if( is.null( g_cent ) == FALSE ){
        g_cent<-as.data.frame( g_cent )
        y2     <- y2 + g_cent[1,1]
        for( i in 1:length(g_cent[-1,1])  ){
          x2[ ,i ]     <- x2[ ,i ] + g_cent[ i + 1, 1 ]
        }

        for( i in 1:length(gc_cent[-1,1])  ){
          xconst2[ ,i ] <- xconst2[ ,i ] + gc_cent[ i + 1, 1 ]
        }
        g_np   <- g_np - 1
      }
        sizelimit2<- sizelimit + g_np + 1
        x2    <- data.frame(x2)
        x2    <- x2[, apply( x2, 2, sd ) > 0 ]
    		res		<- resf_vc_p( y = y2, y0 = y, x = x2, xconst = xconst2, meig = meig, s_id = s_id,
    							penalty = penalty, maxiter = maxiter, sizelimit = sizelimit2, g_np = g_np )
    		b_par	<- res$b
    		s_par	<- res$s
    		sg_par<- NULL
    		e_stat<- res$e
    		b_vc	<- res$b_vc
    		bse_vc<- res$bse_vc
    		t_vc	<- res$t_vc
    		p_vc	<- res$p_vc
    		pred	<- res$pred
    		resid	<- res$resid
    		vc		<- res$vc
    		r		  <- res$r
    		other	<- res$other 

    		const <- 0
    		bg_par<-NULL
    		if( is.null( s_cent ) == FALSE ){
    		  s_ef <- s_cent[, "y" ] - rowSums( data.frame( s_cent [,-( 1:2 )] * b_vc[ t_id == min(as.vector(t_id)), -1] )) - b_vc[ t_id == min(as.vector(t_id)), 1] # t_id is a factor variable, it won't work with the min() routine, so I have to change it to a vector using as.vector() routine. This has to be done with the s_id
    		  if( is.null( sc_cent) == FALSE ){
    		    s_ef<- s_ef - as.matrix( sc_cent [,-( 1:2 )]) %*% b_par[ ,1]
    		  }
    		  s_mean<-mean( s_ef )
    		  s_ef <- s_ef - s_mean
    		  const<- const + s_mean
    		  
    		  bg_s<-data.frame(Estimate = s_ef)
    		  row.names(bg_s) <- paste( "s_", s_id0, sep="")
    		  bg_par<-rbind( bg_par, bg_s)
    		  
    		}
    		if( is.null( t_cent ) == FALSE ){
    		  t_ef <- t_cent[ ,"y" ] - rowSums( data.frame( t_cent [,-( 1:2 )] * b_vc[ s_id == min(as.vector(s_id)), -1 ] )) - b_vc[ s_id == min(as.vector(s_id)), 1]
    		  if( is.null( tc_cent ) == FALSE ){
    		    t_ef<- t_ef - as.matrix( tc_cent [,-( 1:2 )]) %*% b_par[ ,1]
    		  }
    		  t_mean<-mean( t_ef )
    		  t_ef <- t_ef - t_mean
    		  const<- const + t_mean
    		  
    		  bg_t<-data.frame(Estimate = t_ef)
    		  row.names(bg_t) <- paste( "t_", t_id0, sep="")
    		  bg_par <- rbind( bg_par, bg_t)
    		}
    		
    		if( is.null( g_cent ) == FALSE ){
    		  g_ef <- g_cent[ 1, 1 ]
    		  for( cc in 2:(dim(b_vc)[2])){
    		      g_ef <- g_ef - mean( g_cent [ cc, 1 ] * b_vc[, cc ])
    		  }
    		  
    		  if( is.null( gc_cent) == FALSE ){
    		    for( ccc in 1:dim(b_par)[ 1 ])
    		    g_ef <- g_ef - mean( gc_cent[ ccc, 1 ] * b_par[ ccc, 1 ])
    		  }
    		} else {
    		  g_ef  <- 0
    		}
    	  const<- const - g_ef
    	  b_vc[ ,1 ]  <- b_vc[ ,1 ] + const
    	  bse_vc[ ,1 ]<- NA
    	  t_vc[ ,1 ]  <- NA
    		p_vc[ ,1 ]  <- NA
    		
    	} else if( pmodel == "random" ){
    	  g_np      <- 0
    		if( effect == "individual" ){
    			xgroup	<- data.frame( s_id )
    			g_np    <- g_np + length( unique(s_id) ) - 1
        	} else if( effect == "time" ){
    			xgroup	<- data.frame( t_id )
    			g_np    <- g_np + length( unique(t_id) ) - 1
    			
      		} else if( effect == "twoways" ){
    		    xgroup	<- data.frame( s_id, t_id )
    		    g_np    <- g_np + length( unique(s_id)) - 1  + length( unique(t_id)) - 1 
    		}
    	  
    	  sizelimit2<- sizelimit + g_np + 1
		    res		<- resf_vc_p( y = y, x = x, xconst = xconst, xgroup = xgroup, meig = meig, s_id = s_id,
    							penalty = penalty, maxiter = maxiter, sizelimit = sizelimit2 )
    		b_par	<- res$b
    		bg_par	<- res$b_g    	
    		s_par	<- res$s
    		sg_par	<- res$s_g
    		e_stat	<- res$e
    		b_vc		<- res$b_vc
    		bse_vc	<- res$bse_vc
    		t_vc		<- res$t_vc
    		p_vc		<- res$p_vc
    		pred		<- res$pred
    		resid	<- res$resid
    		vc		<- res$vc
    		r		<- res$r
    		other	<- res$other 
    }

    return( list( b = b_par, b_eff = bg_par, s = s_par, s_eff = sg_par, e = e_stat, 
    		b_vc = b_vc, bse_vc = bse_vc, t_vc = t_vc, p_vc = p_vc,
    		pred = pred, resid = resid, vc = vc, r = r, other = other ) )
}




