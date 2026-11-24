#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <cmath>
#include <algorithm>
#include <vector>

// Gaussian profile likelihoods for the full and intercept-only models.
// All vectors are expressed in the eigenbasis of the fixed sampling covariance.
struct Fit { double tau, theta, gamma, obj; };
struct Profile {
  int p,k; const double *y,*v,*q; double scale; bool restricted;
  Fit fixed(double t) const {
    double h00=0,h01=0,h11=0,g0=0,g1=0,ld=0;
    for(int i=0;i<p;i++) {
      double inv=1/(v[i]+t), a=q[i], b=q[p+i];
      h00+=a*a*inv; g0+=a*y[i]*inv;
      if(k==2) {h01+=a*b*inv;h11+=b*b*inv;g1+=b*y[i]*inv;}
      ld+=std::log(v[i]+t);
    }
    double b0=g0/h00,b1=0;
    if(k==2) {double det=h00*h11-h01*h01;
      b0=(h11*g0-h01*g1)/det;b1=(h00*g1-h01*g0)/det;}
    double obj=ld;
    if (restricted) obj += std::log(h00*h11-h01*h01);
    for(int i=0;i<p;i++) {double r=y[i]-q[i]*b0-q[p+i]*b1;obj+=r*r/(v[i]+t);}
    return {t,b0,b1,obj};
  }
  double f(double x) const {return fixed(x*scale).obj;}
  // Golden-section minimization inside each production search bracket.
  // This replaces only the numerical minimizer, not the model or objective.
  double minimize(double lo,double hi) const {
    const double c=(std::sqrt(5.0)-1)/2;
    double x=hi-c*(hi-lo),z=lo+c*(hi-lo),fx=f(x),fz=f(z);
    for(int j=0;j<200 && hi-lo>1e-9;j++) {
      if(fx<=fz) {hi=z;z=x;fz=fx;x=hi-c*(hi-lo);fx=f(x);}
      else {lo=x;x=z;fx=fz;z=lo+c*(hi-lo);fz=f(z);}
    }
    return fx<=fz?x:z;
  }
  Fit run() {
    // Unweighted residual scale, as in method_ML().
    double h00=0,h01=0,h11=0,g0=0,g1=0;
    for(int i=0;i<p;i++){double a=q[i],b=q[p+i];h00+=a*a;g0+=a*y[i];
      h01+=a*b;h11+=b*b;g1+=b*y[i];}
    double b0=g0/h00,b1=0;
    if(k==2){double det=h00*h11-h01*h01;b0=(h11*g0-h01*g1)/det;b1=(h00*g1-h01*g0)/det;}
    double ss=0,sv=0;
    for(int i=0;i<p;i++){double r=y[i]-q[i]*b0-q[p+i]*b1;ss+=r*r;sv+=v[i];}
    scale=std::max(1e-12,std::max(ss/p,sv/p));
    double upper=10,best=0,score=f(0);
    for(int attempt=0;attempt<8;attempt++) {
      std::vector<double> grid(36),val(36);grid[0]=0;val[0]=f(0);
      best=0;score=val[0];
      for(int j=1;j<36;j++){grid[j]=std::exp(std::log(1e-8)+(j-1)*(std::log(upper)-std::log(1e-8))/34);
        val[j]=f(grid[j]);if(val[j]<score){best=grid[j];score=val[j];}}
      for(int j=1;j<35;j++) if(val[j]<=val[j-1] && val[j]<=val[j+1]) {
        double x=minimize(grid[j-1],grid[j+1]),s=f(x);if(s<score){best=x;score=s;}}
      double x=minimize(0,upper),s=f(x);if(s<score){best=x;score=s;}
      if(best<upper*.95) {
        if(f(0)<=f(best)+1e-10)best=0;
        return fixed(best*scale);
      }
      upper*=10;
    }
    Rf_error("Variance profile search exhausted.");
    return fixed(0);
  }
};

extern "C" SEXP lr_profiles(SEXP ys,SEXP vs,SEXP qs) {
  int p=Rf_nrows(ys),n=Rf_ncols(ys);
  if (TYPEOF(ys)!=REALSXP || TYPEOF(vs)!=REALSXP || TYPEOF(qs)!=REALSXP ||
      p<3 || n<1 || Rf_length(vs)!=p || Rf_nrows(qs)!=p || Rf_ncols(qs)!=2)
    Rf_error("Invalid likelihood-ratio inputs.");
  SEXP out=PROTECT(Rf_allocMatrix(REALSXP,n,7));
  for(int b=0;b<n;b++) {
    if(b%100==0) R_CheckUserInterrupt();
    Profile pr{p,1,REAL(ys)+p*b,REAL(vs),REAL(qs),1,false};
    Fit null=pr.run();pr.k=2;Fit full=pr.run();
    double stat=null.obj-full.obj;
    if(!std::isfinite(stat)||stat< -1e-7) Rf_error("Invalid likelihood-ratio statistic.");
    double vals[7]={std::max(0.0,stat),null.tau,full.tau,null.theta,full.theta,full.gamma,full.obj};
    for(int j=0;j<7;j++)REAL(out)[b+n*j]=vals[j];
  }
  UNPROTECT(1);return out;
}

extern "C" SEXP reml_full(SEXP ys,SEXP vs,SEXP qs) {
  int p=Rf_length(ys);
  if(TYPEOF(ys)!=REALSXP || TYPEOF(vs)!=REALSXP || TYPEOF(qs)!=REALSXP ||
     p<3 || Rf_length(vs)!=p || Rf_nrows(qs)!=p || Rf_ncols(qs)!=2)
    Rf_error("Invalid REML inputs.");
  Profile pr{p,2,REAL(ys),REAL(vs),REAL(qs),1,true};
  Fit full=pr.run();
  if(!std::isfinite(full.obj)) Rf_error("Non-finite REML objective.");
  SEXP out=PROTECT(Rf_allocVector(REALSXP,4));
  double vals[4]={full.tau,full.theta,full.gamma,full.obj};
  for(int j=0;j<4;j++) REAL(out)[j]=vals[j];
  UNPROTECT(1); return out;
}
