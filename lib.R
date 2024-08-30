
binned.sample = function(f, goal = 15000, br=30) {
  flims = range(f[-1])
  fc = cut(log(f),breaks = seq(log(flims[1]),log(flims[2]),len=br+1))
  k = ceiling(goal/br)
  i = by(seq_len(nrow(spec)), fc, function(x) if (length(x) < k) x else sample(x, size=k))
  i = sort(do.call(c,i))
}


spec2d = function(tab,span,plot=TRUE,remove_trend=TRUE) {
  taper.w = function(n,p=0.1) {
    w = function(w) {ifelse(abs(w)<1,0.5*(1-cos(w*pi)),1)}
    w(seq_circ(n)/(n*p))
  }
  wx = taper.w(nrow(tab))
  wy = taper.w(ncol(tab))
  tab_taper = tab
  if (remove_trend) {
    data = data.frame(val=as.vector(tab),x=as.vector(row(tab)),y=as.vector(col(tab)))
    m = lm(val~x+y,data=data)
    tab_taper[] = tab_taper - m$fitted.values
  }
  tab_taper = tab_taper * outer(wx,wy)
  tab_p = fft(tab_taper)/length(tab_taper)
  tab_p2 = Re(tab_p*Conj(tab_p))*prod(span)
  tab_f = sqrt(outer((seq_circ(nrow(tab))/span[1])^2,(seq_circ(ncol(tab))/span[2])^2,"+"))
  tab_sp = data.frame( p = as.vector(tab_p2), f = as.vector(tab_f) )
  if (plot) {
    i = binned.sample(tab_sp$f)
    plot(tab_sp$f[i],tab_sp$p[i],log="y")
  }
  invisible(tab_sp)
}
