#include <Rcpp.h>
using namespace Rcpp;

//
// Functions for trajectory analysis
// 
// Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
//
// De Caceres M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (in preparation). Trajectory analysis in community ecology.
//

// Determines the constant that has to be added to each segment
// in order to reach triangle inequality
// [[Rcpp::export(".k2triangleC")]]
double k2triangle(double d1, double d2, double d3) {
  double k = 0.0;
  k = std::max(k, d3 - (d1+d2));
  k = std::max(k, d2 - (d3+d1));
  k = std::max(k, d1 - (d2+d3));
  return(k);
}

//
// Triangle inequality for one triplet
//
// [[Rcpp::export(".triangleinequalityC")]]
bool triangleinequality(double d1, double d2, double d3, double tol=0.0001){
  if((d1+d2)<d3*(1.0-tol)) return(false);
  else if((d1+d3)<d2*(1.0-tol)) return(false);
  else if((d2+d3)<d1*(1.0-tol)) return(false);
  return(true);
}

//
// Determines if the matrix fulfills the triangle inequality
//
// [[Rcpp::export(".ismetricC")]]
bool ismetric(NumericMatrix dmat, double tol=0.0001) {
  int n = dmat.nrow();
  for(int i=0; i<n;i++) {
    for(int j=i;j<n;j++) {
      for(int k=j;k<n; k++) {
        bool ti = triangleinequality(dmat(i,j), dmat(i,k), dmat(j,k), tol);
        if(!ti) return(false);
      }
    }
  }
  return(true);
}

// Projection of one point into a (directed) segment
//
// param dref Distance between the two segment endpoints
// param d1 Distance from the target point to the initial segment endpoint
// param d2 Distance from the target point to the final segment endpoint
//
// returns a vector with:
//   - The distance of between the START of the reference segment and the projected point 
//   - The distance of between the END of the reference segment and the projected point 
//   - The distance height (rejection) of the projection
// [[Rcpp::export(".projectionC")]]
NumericVector projection(double dref, double d1, double d2, bool add = true) {
  if(add) { //Correct triangle inequality if needed
    double k = k2triangle(d1,d2,dref);
    d1 = d1 + k;
    d2 = d2 + k;
    dref = dref + k;
  }
  double a1 = (pow(d1,2.0)+pow(dref,2.0)-pow(d2,2.0))/(2.0*dref);
  double a2 = dref-a1;
  double s = pow(d1,2.0)-pow(a1,2.0);
  double h = NA_REAL;
  if(s>=0.0) h = sqrt(s);
  return(NumericVector::create(a1,a2,h));
}

//
// Angular attribute of a pair of consecutive segments 
// Tripathi, P.K., Debnath, M., & Elmasri, R. 2016. A Direction Based Framework for Trajectory Data Analysis. Proceedings of the 9th ACM International Conference on PErvasive Technologies Related to Assistive Environments - PETRA ’16. 
// doi: 10.1145/2910674.2910728
//
// param d12 Distance from p1 to p2
// param d23 Distance from p2 to p3
// param d13 Distance from p1 to p3
// [[Rcpp::export(".angleConsecutiveC")]]
double angularAttributeConsecutive(double d12, double d23, double d13, bool add = true) {
  if(add) { //Correct triangle inequality if needed
    double k = k2triangle(d12,d23,d13);
    d12 = d12 + k;
    d13 = d13 + k;
    d23 = d23 + k;
  }
  double a1 = (pow(d12,2.0)+pow(d13,2.0)-pow(d23,2.0))/(2.0*d13);
  double a2 = d13-a1;
  if(add) {
    a1 = std::min(a1,d12);
    a2 = std::min(a2,d23);
  }
  double alpha = acos(a1/d12)*(180.0/M_PI);
  double beta = acos(a2/d23)*(180.0/M_PI);
  return(alpha + beta); //Angle between the direction of the first segment and the direction of the second
}

//
// Distance from one point to one segment
//
// [[Rcpp::export(".distanceToSegmentC")]]
NumericVector distanceToSegment(double dref, double d1, double d2, bool add = true) {
  NumericVector p = projection(dref,d1, d2, add);
  if(NumericVector::is_na(p[2]) || (p[0]<0.0) || (p[1]<0.0)) {
    if(d1<d2) {
      p[0] = 0.0;
      p[1] = dref;
      p[2] = d1;
    } else {
      p[0] = dref;
      p[1] = 0.0;
      p[2] = d2;
    }
  }
  return(p);
}

// Based on the law of cosines (https://en.wikipedia.org/wiki/Law_of_cosines)
// For outer triangle we have
// c1 = dref, b1 = d1, a1 = d2
// b1^2 = a1^2 + c1^2 - 2*a1*c1*cos(beta1)
// so that 
//  -2*a1*cos(beta1) = (b1^2 - a1^2 - c1^2/c1
//
// For the inner triangle we have:
// c2 = dref*(1-p), b2 = ?, a2 = d2 
// b2^2 = a2^2 + c2^2 - 2*a2*c2*cos(beta2)
//
// Since a1 = a2 and beta1 = beta2, we have that
//  -2*a1*cos(beta1) = -2*a2*cos(beta2)
// And
// b2^2 = a2^2 + c2*2 - c2*((b1^2 - a1^2 - c1^2/c1)
// [[Rcpp::export(".distanceToInterpolatedC")]]
double distanceToInterpolated(double dref, double d1, double d2, double p, bool add = true) {
  double c2 = dref*(1.0 - p);
  double dsq = (d2*d2) + (c2*c2) + c2*((d1*d1) - (d2*d2)  - (dref*dref))/dref;
  if(add) {
    dsq = std::max(0.0, dsq);
  }
  return(sqrt(dsq));
}

//
// Distance between two segments
//
// [[Rcpp::export(".twoSegmentDistanceC")]]
double twoSegmentDistance(NumericMatrix dmat12, String type="directed-segment", bool add = true) {
  double ds1e1 = dmat12(0,1);
  double ds1s2 = dmat12(0,2);
  double ds1e2 = dmat12(0,3);
  double de1s2 = dmat12(1,2);
  double de1e2 = dmat12(1,3);
  double ds2e2 = dmat12(2,3);
  double Ds = NA_REAL; 
  if((type=="Hausdorff") || (type == "directed-segment")) {
    NumericVector ps1_2 = distanceToSegment(ds2e2,ds1s2, ds1e2, add);
    NumericVector pe1_2 = distanceToSegment(ds2e2,de1s2, de1e2, add);
    NumericVector ps2_1 = distanceToSegment(ds1e1,ds1s2, de1s2, add);
    NumericVector pe2_1 = distanceToSegment(ds1e1,ds1e2, de1e2, add);
    double ds1_2 = ps1_2[2];
    double de1_2 = pe1_2[2];
    double ds2_1 = ps2_1[2];
    double de2_1 = pe2_1[2];
    //Modifications for directionality of segments
    if(type == "directed-segment") { 
      if(ps1_2[0]>pe1_2[0]) {
        de1_2 = std::min(ds1e1+ps1_2[2], ds1e1+pe1_2[2]);
      }
      if(ps2_1[0]>pe2_1[0]) {
        de2_1 = std::min(ds2e2+ps2_1[2],ds2e2+pe2_1[2]);
      }
    }
    Ds = max(NumericVector::create(ds1_2, de1_2, ds2_1, de2_1));
  } else if (type=="PPA"){ //Perpendicular/Parallel/Angle
    if(ds1e1 > ds2e2) { // Switch roles if longest segment is 1
      ds2e2 = dmat12(0,1);
      ds1s2 = dmat12(0,2);
      de1s2 = dmat12(0,3);
      ds1e2 = dmat12(1,2);
      de1e2 = dmat12(1,3);
      ds1e1 = dmat12(2,3);
    }
    //Assumes longer segment is 2
    NumericVector ps1_2 = distanceToSegment(ds2e2,ds1s2, ds1e2, add);
    NumericVector pe1_2 = distanceToSegment(ds2e2,de1s2, de1e2, add);
    double lp1 = ps1_2[2];
    double lp2 = pe1_2[2];
    double dperpendicular = (pow(lp1,2.0)+pow(lp2,2.0))/(lp1+lp2);
    double lpar1 = std::min(ps1_2[0],ps1_2[1]);
    double lpar2 = std::min(pe1_2[0],pe1_2[1]);
    double dparallel = std::min(lpar1,lpar2);
    double dangle = (std::max(lp2,lp1)-std::min(lp2,lp1));
    if(ps1_2[0]>pe1_2[0]) dangle = ds1e1;
          
    Ds = (dperpendicular+dparallel+dangle);
  } 
  return(Ds);
}


//
// Distances between points and arbitrary fuzzy clusters 
//
// param dmat distance matrix (objects in rows and columns)
// param umat membership matrix (objects in rows, clusters in columns)
// [[Rcpp::export(".distanceToClusters")]]
NumericMatrix distanceToClusters(NumericMatrix dmat, NumericMatrix umat) {
  int N = dmat.nrow();
  int K = umat.ncol();
  NumericMatrix d2c(N,K);
  for(int k=0;k<K;k++){
    double cardinality = sum(umat(_,k));
    double vg = 0.0;
    for(int i1=0;i1<N;i1++) {
      for(int i2=0;i2<N;i2++) {
        vg += umat(i1,k)*umat(i2,k)*pow(dmat(i1,i2),2.0);
      }
    }
    vg = vg/(2.0*pow(cardinality, 2.0));
    for(int i1=0;i1<N;i1++) {
      double sqd2c_i1 = 0.0;
      for(int i2=0;i2<N;i2++) {
        sqd2c_i1 += umat(i2,k)*pow(dmat(i1,i2),2.0);
      }
      d2c(i1,k) = sqrt(sqd2c_i1/cardinality - vg);
    }  
  }
  return(d2c);
}

//
// Distances between arbitrary (fuzzy) clusters 
//
// param dmat distance matrix (objects in rows and columns)
// param umat membership matrix (objects in rows, clusters in columns)
// [[Rcpp::export(".distanceBetweenClusters")]]
NumericMatrix distanceBetweenClusters(NumericMatrix dmat, NumericMatrix umat) {
  int N = dmat.nrow();
  int K = umat.ncol();
  NumericVector card(K, 0.0);
  NumericVector var(K, 0.0);
  for(int k=0;k<K;k++){
    card[k] = sum(umat(_,k));
    double vg = 0.0;
    for(int i1=0;i1<N;i1++) {
      for(int i2=0;i2<N;i2++) {
        vg += umat(i1,k)*umat(i2,k)*pow(dmat(i1,i2),2.0);
      }
    }
    var[k] = vg/(2.0*pow(card[k], 2.0));
  }
  NumericMatrix dbc(K,K);
  for(int k1=0;k1<K;k1++){
    dbc(k1, k1) = 0.0;
    for(int k2=0;k2<k1;k2++){
      double cd = 0.0;
      for(int i1=0;i1<N;i1++) {
        for(int i2=0;i2<N;i2++) {
          cd += umat(i1,k1)*pow(dmat(i1,i2),2.0)*umat(i2,k2);
        }
      }
      dbc(k1,k2) = sqrt((cd/(card[k1]*card[k2])) - var[k1] - var[k2]);
      dbc(k2,k1) = dbc(k1,k2);
    }  
  }
  return(dbc);
}

// NOT PRESENTLY USED
double pt(double dIT, double dXT, double dPX, double dPI) {
  if(dPI==0) return(dIT);
  else if(dPX==0) return(dXT);
  double A = pow(dXT,2.0)- pow(dPX,2.0);
  double B = pow(dIT,2.0)- pow(dPI,2.0);
    
  double ax =pow(dXT,2.0)+pow(dIT,2.0);
  double bx = (pow(dXT,2.0)*2.0*B) + (pow(dIT,2.0)*2.0*A) - (4.0*pow(dIT,2.0)*pow(dXT,2.0));
  double cx = pow(dXT,2.0)*pow(B,2.0) + pow(dIT,2.0)*pow(A,2.0);
  double z = pow(bx,2.0)-(4.0*ax*cx);
  double d2 = ((-1.0)*bx + sqrt(z))/(2.0*ax);
  return(sqrt(d2));
}

