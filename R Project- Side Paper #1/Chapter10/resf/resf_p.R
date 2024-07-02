resf_p  	<- function( y, y0 = NULL, x = NULL, xgroup = NULL, s_id = NULL, meig, g_np = NULL, par00 = NULL ){

    lik_resf	<- function( par0, ev, M, m, yy, n, nx, nxx, ne, xg_id ){
    	par	<- par0 ^ 2
    	evv	<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    	evSqrt	<- par[ 2 ] * sqrt( evv )
    	if( is.null( xg_id ) == FALSE ){
    	  for( k in 1:max( xg_id )){
    	    evSqrt<-c( evSqrt, rep( par[ 2 + k ], sum( xg_id == k ) ) ) 
    	  }
    	  neg  <- ne + length( xg_id )
    	} else {
    	  neg   <-ne
    	}
    	Mw	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    	M[  -( 1:nx ), -( 1:nx ) ]	<- Mw + diag( neg )
    	M[  -( 1:nx ),    1:nx   ]	<- M[ -( 1:nx ), 1:nx ] * evSqrt
    	M[     1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ), 1:nx ] )
    	M0				<- M
    	M0[ -( 1:nx ), -( 1:nx ) ]	<- Mw
    	m[  -( 1:nx) ]			<- m[ -( 1:nx ) ] * evSqrt

    	test    <-try( Minv	<- solve( M, tol = 1e-25 ) )
    	if("try-error" %in% class( test )){
    		loglik  <- Inf
    	} else {
    		b	<- Minv %*% m
    		sse	<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
    		dd	<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
    		term1	<- determinant( M )$modulus
    		term2	<- ( n - nxx ) * ( 1 + log( 2 * pi * dd / ( n - nxx ) ) )
    		loglik	<- term1 + term2
    	}
    	return( loglik )
    }

    n		<- length( y )
    if( is.null( x ) ){
    	X	<- as.matrix( rep( 1, n ) )
    	xname	<- "(Intercept)"
    	x_id	<- NULL
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == FALSE ){
    		mode( X00 ) <- "numeric"
    	}
    	x_id	<- apply( X00, 2, sd ) != 0
    	if( sum( x_id ) == 0 ){
    		X	<- as.matrix( rep( 1, n ) )
    		xname	<- "(Intercept)"
    		x_id	<- NULL
    	} else {
    		X0	<- X00[ , x_id ]
    		X	<- as.matrix( cbind( 1, X0 ) )
    		xname	<- c( "(Intercept)", names( as.data.frame( X0 ) ) )
    	}
    }

    if( is.null( xgroup ) == FALSE ){
      xgroup  <- data.frame(xgroup)
      ng	    <- dim( xgroup )[ 2 ]
      xg_id0	<- 1
      for( ff in 1:ng ){
        xg00	  <- factor( xgroup[ , ff ] )
        xg0	    <- model.matrix( ~ 0 + xg00 )
        xg0s	  <- colSums( xg0 )
        xg_idid	<- max( which( xg0s == max( xg0s ) ) )
        Xg0	    <- data.frame( xg0[ , -xg_idid ] )
        names( Xg0 ) <- paste( names( as.data.frame(xgroup) )[ ff ], "_", levels( xg00 )[ -xg_idid ], sep = "" )
        xg_id1	<- rep( xg_id0, length( names( Xg0 ) ) )
        if( ff == 1 ){
          Xg	  <- Xg0
          xg_id	<- xg_id1
        } else {
          Xg	  <- cbind( Xg, Xg0 )
          xg_id	<- c( xg_id, xg_id1 )
        }
        xg_id0	<- xg_id0 + 1
      }
      
    } else {
      ng	<- 0
      xg_id<-NULL
    }
    
    ev		<- meig$ev
    nx		<- dim( X )[ 2 ]
    ne		<- length( ev )
    yy     	<- sum( y ^ 2 )
    XX		<- crossprod(  X )
    Xy		<- crossprod(  X, y )

    if( is.null( g_np ) == FALSE ){
      nxx <- nx + g_np
    } else {
      nxx <- nx
      g_np<- 0
    }
    
    if( is.null( s_id ) == FALSE ){
      meig$sf<-meig$sf[ meig$other$s_id, ][ s_id, ]
      meig$other$fast <- 1
    }
    
    if( ng == 0 ){
      EX		<- crossprod( meig$sf, X )
      Ey		<- crossprod( meig$sf, y )
      if( meig$other$fast == 0 ){
    	  EE	<- diag( ne )
      } else if( meig$other$fast == 1 ){
    	  EE	<- crossprod( meig$sf )
      }
    } else {
      meig$sf <-as.matrix( cbind(meig$sf, Xg) )
      EX		<- crossprod( meig$sf, X )
      Ey		<- crossprod( meig$sf, y )
      EE	  <- crossprod( meig$sf )
    }
    
    M		<- as.matrix( rbind( cbind( XX, t( EX ) ), cbind( EX, EE ) ) )
    m		<- c( Xy, Ey )
    
    if( is.null( par00 )){
      par0<- c( 1, 1 )
      if( ng != 0 ) par0  <- c( par0, rep( 1, ng ) )
      res		<- optim( fn = lik_resf, par0, ev = ev, M = M, m = m, yy = yy,
    	    		   n = n, nx = nx, nxx  = nxx, ne = ne, xg_id = xg_id )
    } else {
      res       <- list(NULL)
      res$par   <- par00
      res$value	<- lik_resf(par0=par00, ev = ev, M = M, m = m, yy = yy,
                      n = n, nx = nx, nxx  = nxx, ne = ne, xg_id = xg_id )
    }
    par		<- res$par ^ 2
    loglik	<- ( -1 / 2 ) * res$value
    evv		<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    evSqrt	<- par[ 2 ] * sqrt( evv )
    
    if( is.null( xg_id ) == FALSE){
      for( k in 1:max( xg_id )){
        evSqrt<-c( evSqrt, rep( par[ 2 + k ], sum( xg_id == k ) ) )
      }
    }
    
    sf2		<- t( t( meig$sf ) * evSqrt )
    X2		<- as.matrix( cbind( X, sf2 ) )
    XX2		<- crossprod( X2 )
    diag( XX2 )[ -( 1:nx ) ] <- diag( XX2 )[ -( 1:nx ) ] + 1
    XXinv2	<- solve( XX2 )
    XXinv2_X	<- XXinv2 %*% t( X2 )
    b		<- XXinv2_X %*% y
    b[ -( 1:nx ) ] <- b[ -( 1:nx ) ] * evSqrt
    pred	<- as.matrix( cbind( X, meig$sf ) ) %*% b
    resid	<- as.numeric(y) - pred # I have to change the y to a numeric vector or it won't calculate because y is a pseries and won't be conformatory with pred, which is a one column matrix.
    SSE		<- sum( resid ^ 2 )
    if( is.null( y0 ) ){
      y0  <- y
    }
    SSY		<- sum( ( y0 - mean( y0 ) ) ^ 2 )
    sig		<- SSE / ( n - nxx )
    bse		<- sqrt( sig ) * sqrt( diag( XXinv2 ) )
    SF		<- meig$sf[, 1:ne ] %*% b[ ( nx + 1 ):( nx + ne ) ]

    np		<- nxx + 1 + 2 + ng
    AIC		<- -2 * loglik + np * 2
    BIC		<- -2 * loglik + np * log( n )
    r2_0	<- 1 - SSE / SSY
    r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1)

    r		  <- b  [ ( nx + 1 ):( nx + ne ) ]
    rse	  <- bse[ ( nx + 1 ):( nx + ne ) ]
    rt	  <- r / rse
    rpar  <- data.frame( Estimate = r, SE = rse, t_value = rt )
    rownames( rpar )	<- paste( "r", 1:ne, sep = "" )
    
    if( ng != 0 ){
      b_g		  <- b  [ -( 1:( nx + ne ) ) ]
      bse_g	  <- bse[ -( 1:( nx + ne ) ) ]
      bt_g	  <- b_g / bse_g
      bpar_g  <- data.frame( Estimate = b_g, SE = bse_g, t_value = bt_g )
    } else {
      bpar_g  <- NULL
    }
    
    bt		<- b[ 1:nx ] / bse[ 1:nx ]
    df		<- sum(t(X2)*XXinv2_X) + g_np
    bp		<- 2 - 2 * pt( abs( bt ), df = n - df )
    b_par	<- data.frame( Estimate = b[ 1:nx ], SE = bse[ 1:nx ], t_value = bt, p_value = bp )
    rownames( b_par ) <- xname

    par[ 2 ]	<- par[ 2 ] * sqrt( sig )
    sf_par	<- data.frame( par = par[ c( 2, 1 ) ] )
    names( sf_par )   <- "Estimate"
    rownames( sf_par )<- c( "shrink_sf_SE", "shrink_sf_alpha" )
    
    if( ng != 0 ){
      parG		<- t( par[ -( 1:2 )] * sqrt( sig ) )
      bg_par	<- data.frame( parG )
      rownames( bg_par )	<- "random_SE"
      names( bg_par )	<- names( as.data.frame( xgroup ) )
    } else {
      bg_par  <-NULL
    }
    
    e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
    
    other	<- list( x_id = x_id, model = "resf", par0 = res$par, nx = nx, df = df )

    return( list( b = b_par, b_g = bpar_g, s = sf_par, s_g = bg_par, e = e_stat,
    		  r = rpar, sf = SF, pred = pred, resid = resid, other = other ) )
}
