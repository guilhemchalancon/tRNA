x <- XPR$expression.noise
require(YET); fNAME( rownames ( x[which( x$dm.sd == max(x$dm.sd,na.rm=T)),] ) )
x <- x[-which( x$dm.sd == max(x$dm.sd,na.rm=T)),]
plot(x$dm.ypd ~ x$dm.sd, pch=16, col="#00000022",xlab="SD medium (DM)", ylab="YPD medium (DM)")
abline(h=mean(x$dm.ypd,na.rm=T),col="orange",lwd=2)
abline(v=mean(x$dm.sd,na.rm=T),col="orange",lwd=2)
abline(0,1,col="red",lwd=2)

# runmed {stats}  R Documentation
# Running Medians – Robust Scatter Plot Smoothing
# Compute running medians of odd span. This is the ‘most robust’ scatter plot smoothing possible. 
# For efficiency (and historical reason), you can use one of two different algorithms giving identical results.

# The two algorithms are internally entirely different:
#   
# "Turlach"
# is the Härdle–Steiger algorithm (see Ref.) as implemented by Berwin Turlach. 
# A tree algorithm is used, ensuring performance O(n * log(k)) where n = length(x) which is asymptotically optimal.
# 
# "Stuetzle"
# is the (older) Stuetzle–Friedman implementation which makes use of median updating when one observation enters and one leaves the smoothing window. 
# While this performs as O(n * k) which is slower asymptotically, it is considerably faster for small k or n.